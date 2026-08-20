#!/usr/bin/env python3

"""Verify tutorial script SHA-256 checksums against the repository manifest."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata" / "script_manifest.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    failures: list[str] = []
    with MANIFEST.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            path = ROOT / row["repository_path"]
            if not path.exists():
                failures.append(f"MISSING: {row['repository_path']}")
                continue
            observed = sha256(path)
            if observed != row["sha256"]:
                failures.append(
                    f"CHECKSUM MISMATCH: {row['repository_path']}\n"
                    f"  expected: {row['sha256']}\n"
                    f"  observed: {observed}"
                )
            else:
                print(f"PASS: {row['repository_path']}")

    if failures:
        raise SystemExit("\n".join(failures))
    print("All script checksums match the manifest.")


if __name__ == "__main__":
    main()