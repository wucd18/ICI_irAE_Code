#!/usr/bin/env python3
"""Resumable multi-range HTTPS downloader for large public GEO files."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import math
import os
import pathlib
import shutil
import time
import urllib.request


def remote_size(url: str) -> tuple[int, str | None, str | None]:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        size = int(response.headers["Content-Length"])
        return size, response.headers.get("ETag"), response.headers.get("Last-Modified")


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(4 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def download_range(
    url: str,
    part_path: pathlib.Path,
    start: int,
    end: int,
    retries: int = 8,
) -> int:
    expected = end - start + 1
    existing = part_path.stat().st_size if part_path.exists() else 0
    if existing > expected:
        raise RuntimeError(f"Oversized partial file: {part_path}")
    if existing == expected:
        return existing

    current = start + existing
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "Range": f"bytes={current}-{end}",
                    "User-Agent": "WCD-irAE-public-data-downloader/1.0",
                },
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                if response.status != 206:
                    raise RuntimeError(
                        f"Server ignored Range request ({response.status}) for {part_path.name}"
                    )
                with part_path.open("ab") as handle:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                        current += len(block)
            final_size = part_path.stat().st_size
            if final_size != expected:
                raise RuntimeError(
                    f"Partial length mismatch for {part_path.name}: {final_size} != {expected}"
                )
            return final_size
        except Exception:
            if attempt == retries:
                raise
            time.sleep(min(2**attempt, 30))
            existing = part_path.stat().st_size if part_path.exists() else 0
            current = start + existing
    raise AssertionError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("output")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-mb", type=int, default=64)
    args = parser.parse_args()

    output = pathlib.Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    total, etag, last_modified = remote_size(args.url)
    if output.exists() and output.stat().st_size == total:
        digest = sha256_file(output)
        print(json.dumps({"status": "already_complete", "bytes": total, "sha256": digest}))
        return

    chunk_bytes = args.chunk_mb * 1024 * 1024
    n_parts = math.ceil(total / chunk_bytes)
    parts = [output.with_name(f"{output.name}.part.{index:04d}") for index in range(n_parts)]

    # Reuse bytes from a previously interrupted single-stream download as part 0.
    if output.exists() and 0 < output.stat().st_size <= min(total, chunk_bytes) and not parts[0].exists():
        output.replace(parts[0])

    tasks = []
    for index, part in enumerate(parts):
        start = index * chunk_bytes
        end = min(total - 1, (index + 1) * chunk_bytes - 1)
        tasks.append((part, start, end))

    started = time.time()
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(download_range, args.url, part, start, end): (part, start, end)
            for part, start, end in tasks
        }
        for future in concurrent.futures.as_completed(future_map):
            future.result()
            completed += 1
            downloaded = sum(p.stat().st_size for p in parts if p.exists())
            elapsed = max(time.time() - started, 0.001)
            print(
                f"parts {completed}/{n_parts}; {downloaded / total:.1%}; "
                f"{downloaded / elapsed / 1024 / 1024:.2f} MiB/s",
                flush=True,
            )

    temporary = output.with_name(f"{output.name}.assembling")
    hasher = hashlib.sha256()
    with temporary.open("wb") as destination:
        for part in parts:
            with part.open("rb") as source:
                while True:
                    block = source.read(4 * 1024 * 1024)
                    if not block:
                        break
                    destination.write(block)
                    hasher.update(block)
    if temporary.stat().st_size != total:
        raise RuntimeError("Assembled file length does not match remote Content-Length")
    os.replace(temporary, output)
    for part in parts:
        part.unlink()

    result = {
        "status": "downloaded",
        "url": args.url,
        "output": str(output),
        "bytes": total,
        "sha256": hasher.hexdigest(),
        "etag": etag,
        "last_modified": last_modified,
        "elapsed_seconds": round(time.time() - started, 2),
    }
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
