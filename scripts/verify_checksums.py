"""Validate the SHA-256 of canonical-LF release files on any platform."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def verify() -> int:
    lines = (ROOT / "checksums.sha256").read_text(encoding="ascii").splitlines()
    for line in lines:
        expected, filename = line.split("  ", 1)
        path = ROOT / filename
        canonical = path.read_bytes().replace(b"\r\n", b"\n")
        actual = hashlib.sha256(canonical).hexdigest()
        if actual.lower() != expected.lower():
            raise ValueError(f"SHA-256 mismatch: {filename}")
    return len(lines)


if __name__ == "__main__":
    print(f"Validated {verify()} canonical-LF SHA-256 entries")
