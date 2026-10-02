#!/usr/bin/env python3
"""Plot per-molecule UMI amplification/read-support distributions.

Celatlas count_detail rows are already collapsed molecule observations after
barcode/UMI correction and gene assignment. The `count` column is the number of
reads supporting each collapsed (Barcode, geneID, UMI) observation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_inputs(count_detail: Path, counts: Path | None, called_only: bool) -> pd.DataFrame:
    data = pd.read_csv(count_detail, sep="\t", usecols=["Barcode", "geneID", "UMI", "count"])
    data["Barcode"] = data["Barcode"].astype(str)
    data["geneID"] = data["geneID"].astype(str)
    data["UMI"] = data["UMI"].astype(str)
    data["count"] = pd.to_numeric(data["count"], errors="coerce").fillna(0).astype(np.int64)
    data = data[data["count"] > 0].copy()

    if called_only:
        if counts is None:
            raise ValueError("--counts is required with --called-only")
        marked = pd.read_csv(counts, sep="\t", usecols=["Barcode", "mark"])
        called = set(marked.loc[marked["mark"] == "CB", "Barcode"].astype(str))
        data = data[data["Barcode"].isin(called)].copy()

    return data


def distribution_table(data: pd.DataFrame) -> pd.DataFrame:
    dist = (
        data.groupby("count", sort=True)
        .size()
        .rename("n_molecules")
        .reset_index()
        .rename(columns={"count": "read_support_per_molecule"})
    )
    dist["total_reads"] = dist["read_support_per_molecule"] * dist["n_molecules"]
    dist["duplicate_reads"] = (dist["read_support_per_molecule"] - 1).clip(lower=0) * dist["n_molecules"]
    total_molecules = dist["n_molecules"].sum()
    total_reads = dist["total_reads"].sum()
    total_duplicates = dist["duplicate_reads"].sum()
    dist["molecule_fraction"] = dist["n_molecules"] / total_molecules if total_molecules else 0.0
    dist["read_fraction"] = dist["total_reads"] / total_reads if total_reads else 0.0
    dist["duplicate_read_fraction"] = dist["duplicate_reads"] / total_duplicates if total_duplicates else 0.0
    return dist


def threshold_summary(counts: np.ndarray) -> pd.DataFrame:
    reads = int(counts.sum())
    molecules = int(len(counts))
    duplicate_reads = int((counts - 1).sum())
    rows = []
    for threshold in [1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000]:
        mask = counts >= threshold
        rows.append(
            {
                "threshold_reads_per_molecule": threshold,
                "n_molecules": int(mask.sum()),
                "molecule_fraction": float(mask.mean()) if molecules else 0.0,
                "read_fraction": float(counts[mask].sum() / reads) if reads else 0.0,
                "duplicate_read_fraction": float((counts[mask] - 1).sum() / duplicate_reads)
                if duplicate_reads
                else 0.0,
            }
        )
    return pd.DataFrame(rows)


def top_contribution(data: pd.DataFrame) -> pd.DataFrame:
    ranked = data.sort_values("count", ascending=False).reset_index(drop=True)
    counts = ranked["count"].to_numpy(dtype=np.int64)
    duplicate = counts - 1
    reads = counts.sum()
    duplicate_reads = duplicate.sum()
    molecules = len(ranked)
    rows = []
    for top_n in [10, 100, 1000, 10000, max(1, int(0.001 * molecules)), max(1, int(0.01 * molecules))]:
        top_n = min(top_n, molecules)
        rows.append(
            {
                "top_n": top_n,
                "min_read_support_in_top_n": int(counts[top_n - 1]),
                "max_read_support": int(counts[0]),
                "read_fraction": float(counts[:top_n].sum() / reads) if reads else 0.0,
                "duplicate_read_fraction": float(duplicate[:top_n].sum() / duplicate_reads)
                if duplicate_reads
                else 0.0,
            }
        )
    return pd.DataFrame(rows)


def quantile_summary(counts: np.ndarray) -> pd.DataFrame:
    rows = []
    for q in [0, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999, 1.0]:
        rows.append({"quantile": q, "read_support_per_molecule": float(np.quantile(counts, q))})
    return pd.DataFrame(rows)


def gene_summary(data: pd.DataFrame) -> pd.DataFrame:
    grouped = data.groupby("geneID").agg(
        n_molecules=("count", "size"),
        total_reads=("count", "sum"),
        median_read_support=("count", "median"),
        max_read_support=("count", "max"),
    )
    grouped["duplicate_reads"] = grouped["total_reads"] - grouped["n_molecules"]
    grouped["duplicate_fraction"] = grouped["duplicate_reads"] / grouped["total_reads"]
    return grouped.sort_values(["max_read_support", "total_reads"], ascending=False).reset_index()


def plot_distribution(dist: pd.DataFrame, output: Path, title: str) -> None:
    fig, ax1 = plt.subplots(figsize=(8, 5))
    x = dist["read_support_per_molecule"].to_numpy()
    y = dist["n_molecules"].to_numpy()
    ax1.bar(x, y, width=0.9, color="#4c78a8", alpha=0.85)
    ax1.set_yscale("log")
    ax1.set_xlabel("Reads supporting one collapsed UMI molecule")
    ax1.set_ylabel("Number of molecules (log scale)")
    ax1.set_title(title)
    ax1.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def plot_cumulative(data: pd.DataFrame, output: Path, title: str) -> None:
    ranked = data.sort_values("count", ascending=False).reset_index(drop=True)
    counts = ranked["count"].to_numpy(dtype=np.float64)
    dup = counts - 1
    x = (np.arange(len(counts)) + 1) / len(counts) * 100.0
    read_cum = np.cumsum(counts) / counts.sum() * 100.0
    dup_cum = np.cumsum(dup) / dup.sum() * 100.0 if dup.sum() else np.zeros_like(dup)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, read_cum, label="Cumulative reads", linewidth=2)
    ax.plot(x, dup_cum, label="Cumulative duplicate reads", linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel("Top molecules ranked by read support (%)")
    ax.set_ylabel("Cumulative contribution (%)")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-detail", required=True)
    parser.add_argument("--counts")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--called-only", action="store_true")
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--no-plots", action="store_true", help="Only write CSV/JSON summaries; skip PDF plots.")
    parser.add_argument(
        "--distribution-only",
        action="store_true",
        help="Write only the read-support distribution PDF, skip the cumulative contribution plot.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    data = read_inputs(
        Path(args.count_detail),
        Path(args.counts) if args.counts else None,
        args.called_only,
    )
    counts = data["count"].to_numpy(dtype=np.int64)
    reads = int(counts.sum())
    molecules = int(len(counts))
    duplicate_reads = int((counts - 1).sum())
    saturation = duplicate_reads / reads * 100.0 if reads else 0.0

    dist = distribution_table(data)
    thresholds = threshold_summary(counts)
    top = top_contribution(data)
    quantiles = quantile_summary(counts)
    genes = gene_summary(data)
    top_rows = data.sort_values("count", ascending=False).head(args.top_n)

    dist.to_csv(outdir / "umi_amplification_distribution.csv", index=False)
    thresholds.to_csv(outdir / "umi_amplification_threshold_summary.csv", index=False)
    top.to_csv(outdir / "umi_amplification_top_contribution.csv", index=False)
    quantiles.to_csv(outdir / "umi_amplification_quantiles.csv", index=False)
    genes.head(200).to_csv(outdir / "umi_amplification_top_genes.csv", index=False)
    top_rows.to_csv(outdir / "top_amplified_umis.csv", index=False)

    summary = {
        "sample": args.sample,
        "called_only": bool(args.called_only),
        "molecules": molecules,
        "reads": reads,
        "duplicate_reads": duplicate_reads,
        "duplicate_saturation_percent": saturation,
        "max_read_support_per_molecule": int(counts.max()) if molecules else 0,
        "median_read_support_per_molecule": float(np.median(counts)) if molecules else 0.0,
        "p99_read_support_per_molecule": float(np.quantile(counts, 0.99)) if molecules else 0.0,
        "p999_read_support_per_molecule": float(np.quantile(counts, 0.999)) if molecules else 0.0,
        "molecules_ge_50_reads": int((counts >= 50).sum()),
        "reads_fraction_from_molecules_ge_50": float(counts[counts >= 50].sum() / reads) if reads else 0.0,
        "duplicate_fraction_from_molecules_ge_50": float(((counts[counts >= 50] - 1).sum()) / duplicate_reads)
        if duplicate_reads
        else 0.0,
    }
    (outdir / "umi_amplification_summary.json").write_text(json.dumps(summary, indent=2))

    if not args.no_plots:
        plot_distribution(
            dist,
            outdir / "umi_amplification_distribution.pdf",
            f"{args.sample}: UMI read-support distribution",
        )
        if not args.distribution_only:
            plot_cumulative(
                data,
                outdir / "umi_amplification_cumulative_contribution.pdf",
                f"{args.sample}: top UMI contribution",
            )

    print(json.dumps(summary, indent=2))
    print(f"Wrote UMI amplification QC to {outdir}")


if __name__ == "__main__":
    main()
