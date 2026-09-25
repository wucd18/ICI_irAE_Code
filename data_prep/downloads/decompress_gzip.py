#!/usr/bin/env python3
"""Stream-decompress a gzip file atomically without loading it into memory."""

from __future__ import annotations

import argparse
import gzip
import os
import pathlib
import shutil


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("destination")
    args = parser.parse_args()

    source = pathlib.Path(args.source).resolve()
    destination = pathlib.Path(args.destination).resolve()
    temporary = destination.with_name(destination.name + ".decompressing")
    if destination.exists() and destination.stat().st_size > 0:
        print(f"already exists: {destination} ({destination.stat().st_size} bytes)")
        return
    with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=8 * 1024 * 1024)
    os.replace(temporary, destination)
    print(f"decompressed: {source.name} -> {destination.name}; {destination.stat().st_size} bytes")


if __name__ == "__main__":
    main()

