#!/usr/bin/env python3
"""Create a small, paired FASTQ demo from a read pair.

The input is read record by record, so the requested number refers to read
pairs rather than compressed bytes.  Both mates are validated and their
headers are replaced with synthetic identifiers; instrument, flow-cell,
lane, tile, coordinates, and sample-index fields are never copied. By default
the biological bases and quality values are replaced with synthetic content.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import TextIO


__version__ = "1.0"


def _open_text(path: Path, mode: str) -> TextIO:
    if str(path).endswith((".gz", ".gzip")):
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def _record(handle: TextIO, path: Path) -> tuple[str, str, str, str] | None:
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        if any(lines[1:]):
            raise ValueError(f"truncated FASTQ record in {path}")
        return None
    if any(not line for line in lines):
        raise ValueError(f"truncated FASTQ record in {path}")
    header, sequence, plus, quality = (line.rstrip("\r\n") for line in lines)
    if not header.startswith("@") or not plus.startswith("+"):
        raise ValueError(f"invalid FASTQ record in {path}: expected @ and + lines")
    if len(sequence) != len(quality):
        raise ValueError(f"sequence/quality length mismatch in {path}")
    return header, sequence, plus, quality


def _write_record(handle: TextIO, mate: int, number: int, sequence: str, quality: str) -> None:
    handle.write(f"@CELATLAS_DEMO:{number:08d}/{mate}\n{sequence}\n+\n{quality}\n")


def _synthetic_sequence(length: int, number: int, mate: int) -> str:
    """Generate deterministic, non-biological bases without copying input."""

    alphabet = b"ACGT"
    output = bytearray()
    block = 0
    while len(output) < length:
        digest = hashlib.sha256(f"CELATLAS_DEMO:{number}:{mate}:{block}".encode("ascii")).digest()
        output.extend(alphabet[value & 3] for value in digest)
        block += 1
    return bytes(output[:length]).decode("ascii")


def _pair_key(header: str) -> str:
    """Return the part of a read name that should be identical for both mates."""

    token = header[1:].split()[0]
    if token.endswith(("/1", "/2")):
        token = token[:-2]
    return token


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def create_demo(r1: Path, r2: Path, outdir: Path, sample: str, reads: int, overwrite: bool, sequence_mode: str = "synthetic") -> dict:
    if reads <= 0:
        raise ValueError("--reads must be a positive integer")
    if sequence_mode not in {"synthetic", "copy"}:
        raise ValueError("--sequence-mode must be synthetic or copy")
    if not r1.is_file() or not r2.is_file():
        raise FileNotFoundError("both --r1 and --r2 must be readable files")
    if r1.resolve() == r2.resolve():
        raise ValueError("--r1 and --r2 must be different files")
    safe_sample = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in sample).strip("._")
    if not safe_sample:
        raise ValueError("--sample must contain at least one filename-safe character")
    outdir.mkdir(parents=True, exist_ok=True)
    out_r1 = outdir / f"{safe_sample}_1.fq.gz"
    out_r2 = outdir / f"{safe_sample}_2.fq.gz"
    manifest_path = outdir / f"{safe_sample}.manifest.json"
    outputs = (out_r1, out_r2, manifest_path)
    if not overwrite and any(path.exists() for path in outputs):
        raise FileExistsError("output exists; use --overwrite to replace it")

    count = 0
    try:
        with _open_text(r1, "r") as in1, _open_text(r2, "r") as in2, gzip.open(out_r1, "wt", encoding="utf-8", newline="") as out1, gzip.open(out_r2, "wt", encoding="utf-8", newline="") as out2:
            for count in range(1, reads + 1):
                rec1 = _record(in1, r1)
                rec2 = _record(in2, r2)
                if rec1 is None or rec2 is None:
                    if rec1 is None and rec2 is None:
                        count -= 1
                        break
                    raise ValueError("R1 and R2 contain different numbers of records")
                header1, seq1, _plus1, qual1 = rec1
                header2, seq2, _plus2, qual2 = rec2
                if _pair_key(header1) != _pair_key(header2):
                    raise ValueError(f"R1/R2 read-name mismatch at pair {count}")
                if sequence_mode == "synthetic":
                    seq1, qual1 = _synthetic_sequence(len(seq1), count, 1), "I" * len(seq1)
                    seq2, qual2 = _synthetic_sequence(len(seq2), count, 2), "I" * len(seq2)
                _write_record(out1, 1, count, seq1, qual1)
                _write_record(out2, 2, count, seq2, qual2)
            else:
                # Detect a short mate mismatch even when --reads is reached.
                extra1, extra2 = _record(in1, r1), _record(in2, r2)
                if (extra1 is None) != (extra2 is None):
                    raise ValueError("R1 and R2 contain different numbers of records")
    except Exception:
        for path in outputs:
            path.unlink(missing_ok=True)
        raise

    result = {
        "format": "Celatlas FASTQ demo manifest",
        "tool": "scripts/create_fastq_demo.py",
        "tool_version": __version__,
        "sample": safe_sample,
        "read_pairs": count,
        "read_length": "variable; copied from input",
        "paired": True,
        "headers_deidentified": True,
        "sequence_mode": sequence_mode,
        "biological_sequences_retained": sequence_mode == "copy",
        "source_paths_recorded": False,
        "files": {
            "r1": {"name": out_r1.name, "sha256": _sha256(out_r1), "bytes": out_r1.stat().st_size},
            "r2": {"name": out_r2.name, "sha256": _sha256(out_r2), "bytes": out_r2.stat().st_size},
        },
    }
    manifest_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True, type=Path, help="Input R1 FASTQ(.gz).")
    parser.add_argument("--r2", required=True, type=Path, help="Input R2 FASTQ(.gz).")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--sample", default="CELATLAS_DEMO", help="Synthetic output sample name.")
    parser.add_argument("--reads", type=int, default=20_000, help="Number of paired reads to copy (default: 20000).")
    parser.add_argument("--sequence-mode", choices=("synthetic", "copy"), default="synthetic", help="synthetic replaces biological sequences and quality values (default); copy retains them and requires data-release authorization.")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output files.")
    args = parser.parse_args(argv)
    try:
        result = create_demo(args.r1, args.r2, args.outdir, args.sample, args.reads, args.overwrite, args.sequence_mode)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
