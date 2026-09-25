#!/usr/bin/env python3
"""Extract selected members from an uncompressed remote TAR using HTTPS ranges."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import pathlib
import re
import urllib.request


BLOCK = 512


def request_range(url: str, start: int, end: int) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "Range": f"bytes={start}-{end}",
            "User-Agent": "WCD-irAE-selective-tar-reader/1.0",
        },
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        if response.status != 206:
            raise RuntimeError(f"Remote server ignored range request: {response.status}")
        return response.read()


def parse_tar_header(header: bytes) -> tuple[str, int]:
    if len(header) != BLOCK:
        raise RuntimeError(f"Expected a 512-byte TAR header, got {len(header)}")
    name = header[0:100].split(b"\0", 1)[0].decode("utf-8")
    prefix = header[345:500].split(b"\0", 1)[0].decode("utf-8")
    if prefix:
        name = f"{prefix}/{name}"
    raw_size = header[124:136].split(b"\0", 1)[0].strip() or b"0"
    size = int(raw_size, 8)
    return name, size


def read_filelist(path: pathlib.Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("#Archive/File") != "File":
                continue
            rows.append((row["Name"], int(row["Size"])))
    if not rows:
        raise RuntimeError("No file entries found in GEO filelist")
    return rows


def member_offsets(entries: list[tuple[str, int]]) -> list[tuple[str, int, int]]:
    offset = 0
    result: list[tuple[str, int, int]] = []
    for name, size in entries:
        result.append((name, size, offset))
        offset += BLOCK + math.ceil(size / BLOCK) * BLOCK
    return result


def extract_one(url: str, target: tuple[str, int, int], output_dir: pathlib.Path) -> dict:
    expected_name, expected_size, header_offset = target
    header = request_range(url, header_offset, header_offset + BLOCK - 1)
    observed_name, observed_size = parse_tar_header(header)
    if observed_name != expected_name or observed_size != expected_size:
        raise RuntimeError(
            "TAR/filelist mismatch at offset "
            f"{header_offset}: expected {expected_name} ({expected_size}), "
            f"observed {observed_name} ({observed_size})"
        )
    data_start = header_offset + BLOCK
    payload = request_range(url, data_start, data_start + expected_size - 1)
    if len(payload) != expected_size:
        raise RuntimeError(f"Short payload for {expected_name}: {len(payload)} != {expected_size}")
    output = output_dir / expected_name
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(output)
    return {
        "member": expected_name,
        "bytes": expected_size,
        "tar_header_offset": header_offset,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "output": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tar_url")
    parser.add_argument("filelist")
    parser.add_argument("member_regex")
    parser.add_argument("output_dir")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    filelist = pathlib.Path(args.filelist).resolve()
    output_dir = pathlib.Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pattern = re.compile(args.member_regex)
    indexed = member_offsets(read_filelist(filelist))
    targets = [entry for entry in indexed if pattern.search(entry[0])]
    if not targets:
        raise RuntimeError("Member regex matched no files")

    results: list[dict] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(extract_one, args.tar_url, target, output_dir): target
            for target in targets
        }
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            print(f"extracted {result['member']} ({result['bytes']} bytes)", flush=True)

    results.sort(key=lambda item: item["member"])
    print(json.dumps({"n_extracted": len(results), "files": results}, ensure_ascii=False))


if __name__ == "__main__":
    main()
