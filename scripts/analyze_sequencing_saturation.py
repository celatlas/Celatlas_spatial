#!/usr/bin/env python3
"""Approximate sequencing saturation and depth sufficiency analysis.

This script works from Celatlas count outputs. It does not reconstruct raw
Cell Ranger read-level barcode/UMI correction, so the saturation estimate is
reported as approximate.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


READ_FRACTIONS = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
READS_PER_CELL_TARGETS = [5_000, 10_000, 15_000, 20_000, 30_000, 40_000, 50_000]


def parse_number(value: str) -> Optional[float]:
    if value is None:
        return None
    match = re.search(r"[-+]?\d[\d,]*(?:\.\d+)?", str(value))
    if not match:
        return None
    return float(match.group(0).replace(",", ""))


def parse_step_metrics(path: Path) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    if not path.exists():
        return metrics
    in_metrics = False
    for line in path.read_text(errors="replace").splitlines():
        if line.strip() == "Metrics":
            in_metrics = True
            continue
        if in_metrics and line.strip() and set(line.strip()) == {"-"}:
            continue
        if in_metrics and line and not line.startswith("- ") and not line.startswith("  "):
            if line.strip().endswith("Metrics") or line.strip() in {"Hidden Metrics", "Help Content"}:
                if line.strip() != "Hidden Metrics":
                    break
        if in_metrics and line.startswith("- ") and ":" in line:
            key, value = line[2:].split(":", 1)
            number = parse_number(value)
            if number is not None:
                metrics[key.strip()] = number
    return metrics


def parse_star_final(path: Path) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    if not path.exists():
        return metrics
    key_map = {
        "Number of input reads": "star_input_reads",
        "Uniquely mapped reads number": "star_uniquely_mapped_reads",
        "Number of reads mapped to multiple loci": "star_multimapped_reads",
        "Number of reads mapped to too many loci": "star_too_many_loci_reads",
    }
    for line in path.read_text(errors="replace").splitlines():
        if "|" not in line:
            continue
        key, value = [part.strip() for part in line.split("|", 1)]
        if key in key_map:
            number = parse_number(value)
            if number is not None:
                metrics[key_map[key]] = number
    mapped = sum(
        metrics.get(key, 0.0)
        for key in ["star_uniquely_mapped_reads", "star_multimapped_reads", "star_too_many_loci_reads"]
    )
    if mapped:
        metrics["star_mapped_reads_unique_multi_toomany"] = mapped
    return metrics


def parse_featurecounts_summary(path: Path) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    if not path.exists():
        return metrics
    df = pd.read_csv(path, sep="\t")
    if df.shape[1] < 2:
        return metrics
    value_col = df.columns[1]
    for _, row in df.iterrows():
        key = str(row["Status"])
        value = parse_number(str(row[value_col]))
        if value is not None:
            metrics[f"featurecounts_{key.lower()}"] = value
    return metrics


def read_count_tables(sample_dir: Path, sample: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    count_detail = sample_dir / "05.count" / f"{sample}_count_detail.txt"
    counts = sample_dir / "05.count" / f"{sample}_counts.txt"
    detail = pd.read_csv(count_detail, sep="\t", usecols=["Barcode", "geneID", "UMI", "count"])
    detail["Barcode"] = detail["Barcode"].astype(str)
    detail["geneID"] = detail["geneID"].astype(str)
    detail["UMI"] = detail["UMI"].astype(str)
    detail["count"] = pd.to_numeric(detail["count"], errors="coerce").fillna(0).astype(np.int64)
    detail = detail[detail["count"] > 0].copy()

    marked = pd.read_csv(counts, sep="\t")
    marked["Barcode"] = marked["Barcode"].astype(str)
    return detail, marked


def summarize_cells(detail_cell: pd.DataFrame, cell_barcodes: Iterable[str]) -> pd.DataFrame:
    cell_index = pd.Index(sorted(cell_barcodes), name="Barcode")
    grouped = detail_cell.groupby("Barcode", sort=False).agg(
        reads=("count", "sum"),
        umi=("UMI", "count"),
        genes=("geneID", "nunique"),
    )
    return grouped.reindex(cell_index, fill_value=0)


def downsample_once(
    data: pd.DataFrame,
    fraction: float,
    rng: np.random.Generator,
    cell_barcodes: Iterable[str],
) -> Dict[str, float]:
    reads = data["count"].to_numpy(dtype=np.int64)
    sampled = rng.binomial(reads, fraction)
    mask = sampled > 0
    sampled_data = data.loc[mask, ["Barcode", "geneID", "UMI"]].copy()
    sampled_data["sampled_count"] = sampled[mask]

    eligible_reads = int(sampled_data["sampled_count"].sum()) if not sampled_data.empty else 0
    dedup_molecules = int(len(sampled_data))
    duplicate_reads = max(0, eligible_reads - dedup_molecules)
    saturation = duplicate_reads / eligible_reads * 100.0 if eligible_reads else 0.0

    cell_summary = summarize_cells(
        sampled_data.rename(columns={"sampled_count": "count"}),
        cell_barcodes,
    )
    cells = len(cell_summary)
    return {
        "eligible_reads": eligible_reads,
        "deduplicated_molecules": dedup_molecules,
        "duplicate_reads": duplicate_reads,
        "duplicate_fraction": saturation / 100.0,
        "sequencing_saturation": saturation,
        "mean_reads_per_cell": eligible_reads / cells if cells else 0.0,
        "median_reads_per_cell": float(cell_summary["reads"].median()) if cells else 0.0,
        "median_umi_per_cell": float(cell_summary["umi"].median()) if cells else 0.0,
        "median_genes_per_cell": float(cell_summary["genes"].median()) if cells else 0.0,
        "total_detected_genes": int(sampled_data["geneID"].nunique()) if not sampled_data.empty else 0,
    }


def downsample_replicates(
    data: pd.DataFrame,
    fractions: List[float],
    repeats: int,
    seed: int,
    cell_barcodes: Iterable[str],
) -> pd.DataFrame:
    records = []
    for fraction in fractions:
        for replicate in range(1, repeats + 1):
            rng = np.random.default_rng(seed + int(round(fraction * 1_000_000)) + replicate)
            result = downsample_once(data, fraction, rng, cell_barcodes)
            result.update(
                {
                    "depth_type": "read_fraction",
                    "read_fraction": fraction,
                    "target_reads_per_cell": np.nan,
                    "replicate": replicate,
                }
            )
            records.append(result)
    return pd.DataFrame(records)


def aggregate_downsampling(replicates: pd.DataFrame) -> pd.DataFrame:
    metric_cols = [
        "eligible_reads",
        "deduplicated_molecules",
        "duplicate_reads",
        "duplicate_fraction",
        "sequencing_saturation",
        "mean_reads_per_cell",
        "median_reads_per_cell",
        "median_umi_per_cell",
        "median_genes_per_cell",
        "total_detected_genes",
    ]
    grouped = replicates.groupby(["depth_type", "read_fraction", "target_reads_per_cell"], dropna=False)
    rows = []
    for keys, part in grouped:
        row = {
            "depth_type": keys[0],
            "read_fraction": keys[1],
            "target_reads_per_cell": keys[2],
            "n_replicates": len(part),
        }
        for col in metric_cols:
            row[f"{col}_mean"] = float(part[col].mean())
            row[f"{col}_sd"] = float(part[col].std(ddof=1)) if len(part) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["depth_type", "read_fraction", "target_reads_per_cell"])


def compute_marginal_gain(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    frac_rows = summary[summary["depth_type"] == "read_fraction"].sort_values("read_fraction")
    prev = None
    for _, row in frac_rows.iterrows():
        if prev is not None:
            delta_reads = row["eligible_reads_mean"] - prev["eligible_reads_mean"]
            delta_umi = row["deduplicated_molecules_mean"] - prev["deduplicated_molecules_mean"]
            delta_median_genes = row["median_genes_per_cell_mean"] - prev["median_genes_per_cell_mean"]
            delta_total_genes = row["total_detected_genes_mean"] - prev["total_detected_genes_mean"]
            rows.append(
                {
                    "from_read_fraction": prev["read_fraction"],
                    "to_read_fraction": row["read_fraction"],
                    "delta_reads": delta_reads,
                    "delta_deduplicated_umi": delta_umi,
                    "umi_gain_per_added_read": delta_umi / delta_reads if delta_reads else np.nan,
                    "delta_median_genes_per_cell": delta_median_genes,
                    "median_gene_gain_per_added_read": delta_median_genes / delta_reads if delta_reads else np.nan,
                    "delta_total_detected_genes": delta_total_genes,
                    "total_gene_gain_per_added_read": delta_total_genes / delta_reads if delta_reads else np.nan,
                    "deduplicated_umi_percent_gain": delta_umi / prev["deduplicated_molecules_mean"] * 100.0
                    if prev["deduplicated_molecules_mean"]
                    else np.nan,
                    "median_genes_percent_gain": delta_median_genes / prev["median_genes_per_cell_mean"] * 100.0
                    if prev["median_genes_per_cell_mean"]
                    else np.nan,
                }
            )
        prev = row
    return pd.DataFrame(rows)


def add_uniform_depth_rows(
    data: pd.DataFrame,
    base_replicates: pd.DataFrame,
    targets: List[int],
    repeats: int,
    seed: int,
    cell_barcodes: Iterable[str],
) -> pd.DataFrame:
    full_reads = int(data["count"].sum())
    n_cells = len(list(cell_barcodes))
    records = []
    for target in targets:
        target_reads = target * n_cells
        if target_reads > full_reads:
            continue
        fraction = target_reads / full_reads
        for replicate in range(1, repeats + 1):
            rng = np.random.default_rng(seed + target + replicate * 17)
            result = downsample_once(data, fraction, rng, cell_barcodes)
            result.update(
                {
                    "depth_type": "target_reads_per_cell",
                    "read_fraction": fraction,
                    "target_reads_per_cell": target,
                    "replicate": replicate,
                }
            )
            records.append(result)
    if records:
        return pd.concat([base_replicates, pd.DataFrame(records)], ignore_index=True)
    return base_replicates


def write_curve(
    summary: pd.DataFrame,
    y_col: str,
    ylabel: str,
    title: str,
    output: Path,
) -> None:
    data = summary[summary["depth_type"] == "read_fraction"].sort_values("mean_reads_per_cell_mean")
    fig, ax = plt.subplots(figsize=(7, 5))
    x = data["mean_reads_per_cell_mean"]
    y = data[y_col]
    yerr_col = y_col.replace("_mean", "_sd")
    yerr = data[yerr_col] if yerr_col in data.columns else None
    ax.errorbar(x, y, yerr=yerr, marker="o", linewidth=1.8, capsize=3)
    ax.set_xlabel("Mean eligible molecule-supporting reads per cell")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def write_summary_text(
    output: Path,
    sample: str,
    summary_row: Dict[str, object],
    marginal: pd.DataFrame,
) -> None:
    last20 = marginal[(marginal["from_read_fraction"] == 0.8) & (marginal["to_read_fraction"] == 0.9)]
    last10 = marginal[(marginal["from_read_fraction"] == 0.9) & (marginal["to_read_fraction"] == 1.0)]
    from80 = marginal[marginal["from_read_fraction"] >= 0.8]
    umi_gain_80_100 = float(from80["delta_deduplicated_umi"].sum()) if not from80.empty else np.nan
    median_gene_gain_80_100 = float(from80["delta_median_genes_per_cell"].sum()) if not from80.empty else np.nan

    saturation = float(summary_row["sequencing_saturation_percent"])
    median_genes = float(summary_row["median_genes_per_cell"])
    median_umi = float(summary_row["median_umi_per_cell"])
    plateau = "yes" if median_gene_gain_80_100 < max(1.0, median_genes * 0.05) else "partial/no"

    text = f"""# Sequencing Saturation And Depth Sufficiency Summary

Sample: `{sample}`

This analysis is marked **approximate sequencing saturation**. It uses Celatlas
`05.count` molecule-supporting read counts after pipeline barcode correction,
UMI correction, multigene UMI resolution, and called-barcode filtering. It does
not replay raw FASTQ read-level Cell Ranger correction/deduplication rules.

## Current Depth

- Approximate sequencing saturation: {saturation:.2f}%
- Eligible molecule-supporting reads in called barcodes: {int(summary_row['eligible_molecule_supporting_reads']):,}
- Deduplicated molecules / UMIs: {int(summary_row['deduplicated_molecules']):,}
- Duplicate reads: {int(summary_row['duplicate_reads']):,}
- Mean eligible reads per cell: {float(summary_row['mean_eligible_reads_per_cell']):.2f}
- Median UMI per cell: {median_umi:.2f}
- Median genes per cell: {median_genes:.2f}

## Depth Sufficiency

From 80% to 100% of current eligible reads, the approximate downsampling curve
adds {umi_gain_80_100:,.0f} deduplicated UMIs and {median_gene_gain_80_100:.2f}
median genes per cell. Platform call: **{plateau}**.

Interpretation: sequencing saturation is a duplicate-read fraction among
eligible reads. It is not the fraction of true RNA molecules captured. Whether
additional sequencing is worthwhile should be judged from the UMI and gene
rarefaction curves and marginal gains.

## Notes

- Raw FASTQ reads are not used as the saturation denominator.
- featureCounts `Assigned` can exceed STAR input reads because it counts
  feature-assignment/alignment records in this pipeline; it is reported as an
  upstream assignment metric, not as the saturation denominator.
- Multi-sample comparison should be repeated at a common reads/cell depth.
"""
    output.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()

    sample_dir = Path(args.sample_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    detail, marked = read_count_tables(sample_dir, args.sample)
    cell_barcodes = set(marked.loc[marked["mark"] == "CB", "Barcode"].astype(str))
    detail_cell = detail[detail["Barcode"].isin(cell_barcodes)].copy()
    cell_summary = summarize_cells(detail_cell, cell_barcodes)

    barcode_metrics = parse_step_metrics(sample_dir / "01.barcode" / "step.log")
    count_metrics = parse_step_metrics(sample_dir / "05.count" / "step.log")
    star_metrics = parse_star_final(sample_dir / "03.star" / f"{args.sample}_Log.final.out")
    featurecounts_metrics = parse_featurecounts_summary(sample_dir / "04.featureCounts" / f"{args.sample}.summary")

    eligible_reads = int(detail_cell["count"].sum())
    dedup_molecules = int(len(detail_cell))
    duplicate_reads = eligible_reads - dedup_molecules
    saturation = duplicate_reads / eligible_reads * 100.0 if eligible_reads else 0.0
    all_transcriptome_reads = int(detail["count"].sum())
    fraction_reads_in_cells = eligible_reads / all_transcriptome_reads * 100.0 if all_transcriptome_reads else 0.0
    n_cells = len(cell_barcodes)

    summary_row = {
        "sample": args.sample,
        "saturation_status": "approximate",
        "saturation_numerator_definition": "eligible_molecule_supporting_reads_in_called_barcodes - deduplicated_molecules",
        "saturation_denominator_definition": "eligible molecule-supporting reads in called barcodes from 05.count count_detail",
        "raw_fastq_reads": int(barcode_metrics.get("Raw Reads", np.nan)) if "Raw Reads" in barcode_metrics else np.nan,
        "valid_barcode_reads": int(barcode_metrics.get("Valid Reads", np.nan)) if "Valid Reads" in barcode_metrics else np.nan,
        "star_input_reads_after_barcode_filtering": int(star_metrics.get("star_input_reads", np.nan))
        if "star_input_reads" in star_metrics
        else np.nan,
        "star_uniquely_mapped_reads": int(star_metrics.get("star_uniquely_mapped_reads", np.nan))
        if "star_uniquely_mapped_reads" in star_metrics
        else np.nan,
        "star_mapped_reads_unique_multi_toomany": int(star_metrics.get("star_mapped_reads_unique_multi_toomany", np.nan))
        if "star_mapped_reads_unique_multi_toomany" in star_metrics
        else np.nan,
        "featurecounts_assigned_records_not_saturation_denominator": int(featurecounts_metrics.get("featurecounts_assigned", np.nan))
        if "featurecounts_assigned" in featurecounts_metrics
        else np.nan,
        "transcriptome_molecule_supporting_reads_all_barcodes": all_transcriptome_reads,
        "confidently_gene_assigned_molecule_supporting_reads_called_barcodes": eligible_reads,
        "eligible_molecule_supporting_reads": eligible_reads,
        "deduplicated_molecules": dedup_molecules,
        "duplicate_reads": duplicate_reads,
        "duplicate_fraction": duplicate_reads / eligible_reads if eligible_reads else 0.0,
        "sequencing_saturation_percent": saturation,
        "estimated_cells": n_cells,
        "mean_eligible_reads_per_cell": eligible_reads / n_cells if n_cells else 0.0,
        "pipeline_valid_reads_per_cell": int(count_metrics.get("Mean Reads per Cell", np.nan))
        if "Mean Reads per Cell" in count_metrics
        else np.nan,
        "median_reads_per_cell": float(cell_summary["reads"].median()) if n_cells else 0.0,
        "median_umi_per_cell": float(cell_summary["umi"].median()) if n_cells else 0.0,
        "median_genes_per_cell": float(cell_summary["genes"].median()) if n_cells else 0.0,
        "total_detected_genes": int(detail_cell["geneID"].nunique()),
        "fraction_reads_in_cells_percent": fraction_reads_in_cells,
    }

    replicates = downsample_replicates(detail_cell, READ_FRACTIONS, args.repeats, args.seed, cell_barcodes)
    replicates = add_uniform_depth_rows(
        detail_cell,
        replicates,
        READS_PER_CELL_TARGETS,
        args.repeats,
        args.seed + 10_000,
        cell_barcodes,
    )
    downsample_summary = aggregate_downsampling(replicates)
    marginal = compute_marginal_gain(downsample_summary)

    pd.DataFrame([summary_row]).to_csv(outdir / "sequencing_saturation_summary.csv", index=False)
    downsample_summary.to_csv(outdir / "sequencing_downsampling_summary.csv", index=False)
    replicates.to_csv(outdir / "sequencing_downsampling_replicates.csv", index=False)
    marginal.to_csv(outdir / "sequencing_marginal_gain.csv", index=False)

    write_curve(
        downsample_summary,
        "sequencing_saturation_mean",
        "Approximate sequencing saturation (%)",
        f"{args.sample}: sequencing saturation vs depth",
        outdir / "sequencing_saturation_curve.pdf",
    )
    write_curve(
        downsample_summary,
        "median_umi_per_cell_mean",
        "Median UMI per cell",
        f"{args.sample}: UMI rarefaction",
        outdir / "umi_rarefaction_curve.pdf",
    )
    write_curve(
        downsample_summary,
        "median_genes_per_cell_mean",
        "Median genes per cell",
        f"{args.sample}: gene rarefaction",
        outdir / "gene_rarefaction_curve.pdf",
    )
    write_curve(
        downsample_summary,
        "total_detected_genes_mean",
        "Total detected genes",
        f"{args.sample}: total detected genes vs depth",
        outdir / "total_detected_genes_curve.pdf",
    )

    params = {
        "sample_dir": str(sample_dir),
        "sample": args.sample,
        "outdir": str(outdir),
        "read_fractions": READ_FRACTIONS,
        "reads_per_cell_targets": READS_PER_CELL_TARGETS,
        "repeats": args.repeats,
        "seed": args.seed,
        "saturation_status": "approximate",
    }
    (outdir / "sequencing_saturation_parameters.json").write_text(json.dumps(params, indent=2))
    write_summary_text(outdir / "sequencing_saturation_interpretation.md", args.sample, summary_row, marginal)

    print(f"Wrote sequencing saturation analysis to {outdir}")


if __name__ == "__main__":
    main()
