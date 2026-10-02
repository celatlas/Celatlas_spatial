#!/usr/bin/env python3
"""Validate that a wheel or source archive contains only public release data."""

from __future__ import annotations

import argparse
import re
import sys
import tarfile
import zipfile
from pathlib import Path


FORBIDDEN_MEMBER_PATTERNS = (
    re.compile(r"(^|/)private(/|$)"),
    re.compile(r"(^|/)research_project_info(?:_|\.)", re.IGNORECASE),
    re.compile(r"(^|/)logs(/|$)"),
    re.compile(r"(^|/)tmp(/|$)"),
    re.compile(r"\.xlsx$", re.IGNORECASE),
)

FORBIDDEN_TEXT = (
    "/home/bioinfo/",
    "/home/smart_prod/",
    "/mnt/service/",
    "research_project_info_new.xlsx",
)


def _members(path: Path) -> list[tuple[str, bytes | None]]:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            return [(info.filename, archive.read(info)) for info in archive.infolist()]
    if path.name.endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(path) as archive:
            result = []
            for info in archive.getmembers():
                payload = archive.extractfile(info).read() if info.isfile() else None
                result.append((info.name, payload))
            return result
    raise ValueError(f"Unsupported release artifact: {path}")


def check_artifact(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        members = _members(path)
    except (OSError, ValueError, zipfile.BadZipFile, tarfile.TarError) as exc:
        return [f"{path}: cannot read archive: {exc}"]

    for name, payload in members:
        if any(pattern.search(name) for pattern in FORBIDDEN_MEMBER_PATTERNS):
            errors.append(f"{path}: forbidden member {name}")
        if name.endswith("scripts/check_public_release.py"):
            continue
        if payload is None or not name.endswith((".py", ".sh", ".txt", ".cfg", ".in", ".yaml", ".yml", ".md")):
            continue
        text = payload.decode("utf-8", errors="replace")
        for marker in FORBIDDEN_TEXT:
            if marker in text:
                errors.append(f"{path}: forbidden text {marker!r} in {name}")
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args(argv)
    errors = [error for artifact in args.artifacts for error in check_artifact(artifact)]
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    for artifact in args.artifacts:
        print(f"Public release audit passed: {artifact}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
