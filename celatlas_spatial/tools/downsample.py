"""
Deterministic downsampling metrics for sequencing saturation plots.
"""

from typing import Iterable, Optional

import numpy as np
import pandas as pd


READ_FRACTION = "read_fraction"
MEDIAN_GENE_NUMBER = "median_gene_number"
UMI_SATURATION = "umi_saturation"
READ_SATURATION = "read_saturation"
ORIGINAL_BARCODE = "OriginalBarcode"

DEFAULT_READ_FRACTIONS = [0.0] + [round(x, 1) for x in np.arange(0.1, 1.1, 0.1)]


def _observed_probability(read_counts: np.ndarray, fraction: float) -> np.ndarray:
    """Expected probability that an observation with read count c is seen at fraction f."""
    if fraction <= 0:
        return np.zeros_like(read_counts, dtype=np.float64)
    if fraction >= 1:
        return np.ones_like(read_counts, dtype=np.float64)

    counts = read_counts.astype(np.float64, copy=False)
    return -np.expm1(counts * np.log1p(-fraction))


def deterministic_downsample_metrics(
    count_detail: pd.DataFrame,
    fractions: Optional[Iterable[float]] = None,
    barcode_col: str = "Barcode",
    gene_col: str = "geneID",
    umi_col: str = "UMI",
    count_col: str = "count",
    molecule_barcode_col: Optional[str] = None,
) -> pd.DataFrame:
    """Calculate deterministic expected downsampling metrics.

    The saturation curve uses the expected number of observed deduplicated
    molecule barcode-gene-UMI observations at each read fraction instead of one
    random subsample. Median genes are calculated from the expected
    gene-detection count per displayed barcode/bin at the same fraction.

    Spatial binning can combine many original capture barcodes into one bin.
    In that case, pass molecule_barcode_col, or include an OriginalBarcode
    column, so identical UMI sequences from different original barcodes are not
    collapsed into one molecule during saturation estimation.
    """
    fractions = list(DEFAULT_READ_FRACTIONS if fractions is None else fractions)
    required_cols = [barcode_col, gene_col, umi_col, count_col]
    if molecule_barcode_col is None and ORIGINAL_BARCODE in count_detail.columns and ORIGINAL_BARCODE != barcode_col:
        molecule_barcode_col = ORIGINAL_BARCODE
    if molecule_barcode_col is not None:
        required_cols.append(molecule_barcode_col)

    required_cols = list(dict.fromkeys(required_cols))
    missing_cols = [col for col in required_cols if col not in count_detail.columns]
    if missing_cols:
        raise ValueError(f"count_detail missing required columns: {', '.join(missing_cols)}")

    data = count_detail[required_cols].copy()
    data[count_col] = pd.to_numeric(data[count_col], errors="coerce").fillna(0)
    data = data[data[count_col] > 0]
    if molecule_barcode_col is not None:
        data[molecule_barcode_col] = data[molecule_barcode_col].fillna(data[barcode_col])

    results = {
        READ_FRACTION: [],
        MEDIAN_GENE_NUMBER: [],
        UMI_SATURATION: [],
        READ_SATURATION: [],
    }
    if data.empty:
        for fraction in fractions:
            results[READ_FRACTION].append(round(float(fraction), 4))
            results[MEDIAN_GENE_NUMBER].append(0.0)
            results[UMI_SATURATION].append(0.0)
            results[READ_SATURATION].append(0.0)
        return pd.DataFrame(results)

    molecule_group_cols = [molecule_barcode_col or barcode_col, gene_col, umi_col]
    umi_counts = (
        data.groupby(molecule_group_cols, sort=False, observed=True)[count_col]
        .sum()
        .to_numpy(dtype=np.float64)
    )
    total_reads = float(umi_counts.sum())

    gene_counts = (
        data.groupby([barcode_col, gene_col], sort=False, observed=True)[count_col]
        .sum()
        .reset_index()
    )
    gene_read_counts = gene_counts[count_col].to_numpy(dtype=np.float64)

    for fraction in fractions:
        fraction = float(fraction)
        if fraction <= 0 or total_reads <= 0:
            saturation = 0.0
            median_genes = 0.0
        else:
            expected_reads = fraction * total_reads
            expected_deduped = float(_observed_probability(umi_counts, fraction).sum())
            saturation = (1.0 - expected_deduped / expected_reads) * 100.0
            saturation = float(np.clip(saturation, 0.0, 100.0))

            gene_counts["_expected_detected"] = _observed_probability(gene_read_counts, fraction)
            expected_genes = gene_counts.groupby(barcode_col, sort=False, observed=True)["_expected_detected"].sum()
            median_genes = float(expected_genes.median()) if not expected_genes.empty else 0.0

        saturation = round(saturation, 2)
        results[READ_FRACTION].append(round(fraction, 4))
        results[MEDIAN_GENE_NUMBER].append(round(median_genes, 2))
        results[UMI_SATURATION].append(saturation)
        results[READ_SATURATION].append(saturation)

    return pd.DataFrame(
        results,
        columns=[READ_FRACTION, MEDIAN_GENE_NUMBER, UMI_SATURATION, READ_SATURATION],
    )
