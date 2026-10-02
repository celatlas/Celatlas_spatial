#!/usr/bin/env python3
"""Convert a GFF3 annotation to a GTF usable by Celatlas RNA steps.

The converter is intentionally conservative: gene, transcript/mRNA, exon, CDS,
and UTR features are preserved, and child features receive both gene_id and
transcript_id attributes. This fixes common plant GFF3 annotations where exons
only carry Parent=transcript_id, which causes featureCounts and downstream
Celatlas GTF parsers to fail when they require gene_id.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote


TRANSCRIPT_FEATURES = {"mRNA", "transcript"}
PASSTHROUGH_CHILD_FEATURES = {
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "UTR",
}


@dataclass(frozen=True)
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: str
    end: str
    score: str
    strand: str
    phase: str
    attrs: dict[str, str]


def normalize_id(value: str | None) -> str:
    if not value:
        return ""
    value = value.strip()
    for prefix in ("gene:", "transcript:", "rna:", "mRNA:"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def parse_attrs(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for field in raw.strip().split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
        elif " " in field:
            key, value = field.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = unquote(value.strip())
    return attrs


def parse_gff3(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) != 9:
                raise ValueError(f"{path}:{line_no}: expected 9 GFF3 columns, got {len(row)}")
            features.append(
                Feature(
                    seqid=row[0],
                    source=row[1],
                    feature_type=row[2],
                    start=row[3],
                    end=row[4],
                    score=row[5],
                    strand=row[6],
                    phase=row[7],
                    attrs=parse_attrs(row[8]),
                )
            )
    return features


def gtf_attrs(attrs: dict[str, str]) -> str:
    fields = []
    for key, value in attrs.items():
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        fields.append(f'{key} "{escaped}";')
    return " ".join(fields)


def parent_ids(attrs: dict[str, str]) -> list[str]:
    raw = attrs.get("Parent", "")
    return [normalize_id(item) for item in raw.split(",") if item.strip()]


def first_value(attrs: dict[str, str], *keys: str, default: str = "") -> str:
    for key in keys:
        value = attrs.get(key)
        if value:
            return normalize_id(value)
    return default


def build_maps(features: Iterable[Feature]) -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    genes: dict[str, dict[str, str]] = {}
    transcript_to_gene: dict[str, str] = {}

    for feature in features:
        if feature.feature_type == "gene":
            gene_id = first_value(feature.attrs, "ID", "Name")
            if not gene_id:
                continue
            genes[gene_id] = {
                "gene_id": gene_id,
                "gene_name": first_value(feature.attrs, "Name", "gene_name", default=gene_id),
                "gene_biotype": first_value(
                    feature.attrs,
                    "biotype",
                    "gene_biotype",
                    "gene_type",
                    default="",
                ),
            }

    for feature in features:
        if feature.feature_type not in TRANSCRIPT_FEATURES:
            continue
        transcript_id = first_value(feature.attrs, "ID", "Name")
        parents = parent_ids(feature.attrs)
        if transcript_id and parents:
            transcript_to_gene[transcript_id] = parents[0]

    return genes, transcript_to_gene


def feature_attrs(
    feature: Feature,
    gene_id: str,
    transcript_id: str = "",
    genes: dict[str, dict[str, str]] | None = None,
) -> dict[str, str]:
    gene_meta = (genes or {}).get(gene_id, {})
    attrs: dict[str, str] = {
        "gene_id": gene_id,
        "gene_name": gene_meta.get("gene_name") or first_value(feature.attrs, "Name", default=gene_id),
    }
    gene_biotype = gene_meta.get("gene_biotype") or first_value(
        feature.attrs,
        "biotype",
        "gene_biotype",
        "gene_type",
        default="",
    )
    if gene_biotype:
        attrs["gene_biotype"] = gene_biotype
    if transcript_id:
        attrs["transcript_id"] = transcript_id
    return attrs


def converted_rows(features: list[Feature]) -> Iterable[list[str]]:
    genes, transcript_to_gene = build_maps(features)
    emitted_gene_ids: set[str] = set()

    for feature in features:
        if feature.feature_type == "gene":
            gene_id = first_value(feature.attrs, "ID", "Name")
            if not gene_id:
                continue
            emitted_gene_ids.add(gene_id)
            attrs = feature_attrs(feature, gene_id, genes=genes)
            yield [
                feature.seqid,
                feature.source,
                "gene",
                feature.start,
                feature.end,
                feature.score,
                feature.strand,
                ".",
                gtf_attrs(attrs),
            ]
            continue

        if feature.feature_type in TRANSCRIPT_FEATURES:
            transcript_id = first_value(feature.attrs, "ID", "Name")
            parents = parent_ids(feature.attrs)
            gene_id = parents[0] if parents else transcript_to_gene.get(transcript_id, "")
            if not gene_id or not transcript_id:
                continue
            attrs = feature_attrs(feature, gene_id, transcript_id, genes)
            yield [
                feature.seqid,
                feature.source,
                "transcript",
                feature.start,
                feature.end,
                feature.score,
                feature.strand,
                ".",
                gtf_attrs(attrs),
            ]
            continue

        if feature.feature_type not in PASSTHROUGH_CHILD_FEATURES:
            continue

        parents = parent_ids(feature.attrs)
        for transcript_id in parents:
            gene_id = transcript_to_gene.get(transcript_id, transcript_id)
            if not gene_id:
                continue
            if gene_id not in emitted_gene_ids and gene_id in genes:
                emitted_gene_ids.add(gene_id)
            attrs = feature_attrs(feature, gene_id, transcript_id, genes)
            yield [
                feature.seqid,
                feature.source,
                feature.feature_type,
                feature.start,
                feature.end,
                feature.score,
                feature.strand,
                feature.phase if feature.feature_type == "CDS" else ".",
                gtf_attrs(attrs),
            ]


def write_gtf(in_gff3: Path, out_gtf: Path) -> None:
    features = parse_gff3(in_gff3)
    out_gtf.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with out_gtf.open("w", newline="") as handle:
        for row in converted_rows(features):
            handle.write("\t".join(row))
            handle.write("\n")
            row_count += 1
    print(f"Wrote {row_count:,} GTF records: {out_gtf}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert GFF3 to a Celatlas-compatible GTF with gene_id on gene/exon/CDS records."
    )
    parser.add_argument("gff3", type=Path, help="Input GFF3 annotation")
    parser.add_argument("gtf", type=Path, help="Output normalized GTF")
    args = parser.parse_args()

    write_gtf(args.gff3, args.gtf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
