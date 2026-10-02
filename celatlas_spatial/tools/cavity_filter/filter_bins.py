"""Filter square-bin 10X matrices with a binary tissue/cavity mask."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy import io as sio

Image.MAX_IMAGE_PIXELS = None


def load_barcodes(barcodes_path: Path) -> list[str]:
    with gzip.open(barcodes_path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def load_features(features_path: Path) -> str:
    with gzip.open(features_path, "rt") as handle:
        return handle.read()


def load_positions(positions_path: Path) -> list[list[str]]:
    rows = []
    with open(positions_path, "r", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if row:
                rows.append(row)
    return rows


def load_mask(mask_path: Path, threshold: int = 127, remove_mask_value: str = "black", dilate: int = 0) -> np.ndarray:
    try:
        mask = np.asarray(Image.open(mask_path).convert("L"), dtype=np.uint8)
    except Exception as exc:
        raise FileNotFoundError(f"Failed to read mask image: {mask_path}")
    binary = mask > int(threshold)
    if remove_mask_value == "black":
        binary = ~binary
    remove_mask = binary.astype(bool)
    if dilate > 0:
        remove_mask = ndi.binary_dilation(remove_mask, iterations=int(dilate))
    return remove_mask


def get_position_extent(position_rows: list[list[str]]) -> tuple[float, float]:
    max_row = max(float(row[4]) for row in position_rows)
    max_col = max(float(row[5]) for row in position_rows)
    return max_row, max_col


def maybe_resize_mask_to_hires(
    remove_mask: np.ndarray,
    hires_image_path: Path,
    position_rows: list[list[str]],
    mask_space: str = "auto",
) -> tuple[np.ndarray, dict | None]:
    if not hires_image_path.exists():
        return remove_mask, None
    try:
        hires = Image.open(hires_image_path)
    except Exception:
        return remove_mask, None

    hires_shape = (hires.size[1], hires.size[0])
    mask_shape = remove_mask.shape[:2]
    max_row, max_col = get_position_extent(position_rows)
    coord_fits_hires = max_row <= hires_shape[0] * 1.05 and max_col <= hires_shape[1] * 1.05
    mask_matches_hires = mask_shape == hires_shape

    should_resize = False
    if mask_space == "hires":
        should_resize = not mask_matches_hires
    elif mask_space == "auto":
        should_resize = coord_fits_hires and not mask_matches_hires

    info = {
        "mask_shape": list(mask_shape),
        "hires_shape": list(hires_shape),
        "resized_to_hires": bool(should_resize),
        "max_position_row": float(max_row),
        "max_position_col": float(max_col),
    }
    if not should_resize:
        return remove_mask, info

    try:
        nearest = Image.Resampling.NEAREST
    except AttributeError:
        nearest = Image.NEAREST
    resized = np.asarray(
        Image.fromarray(remove_mask.astype(np.uint8) * 255).resize((hires_shape[1], hires_shape[0]), nearest)
    ) > 0
    return resized, info


def should_remove(position_row: list[str], remove_mask: np.ndarray, radius: int = 0) -> bool:
    height, width = remove_mask.shape
    px_row = int(round(float(position_row[4])))
    px_col = int(round(float(position_row[5])))
    if px_row < 0 or px_row >= height or px_col < 0 or px_col >= width:
        return False
    if radius <= 0:
        return bool(remove_mask[px_row, px_col])
    y0 = max(0, px_row - radius)
    y1 = min(height, px_row + radius + 1)
    x0 = max(0, px_col - radius)
    x1 = min(width, px_col + radius + 1)
    return bool(remove_mask[y0:y1, x0:x1].any())


def write_barcodes(barcodes: list[str], output_path: Path) -> None:
    with gzip.open(output_path, "wt") as handle:
        for barcode in barcodes:
            handle.write(f"{barcode}\n")


def write_features(features_text: str, output_path: Path) -> None:
    with gzip.open(output_path, "wt") as handle:
        handle.write(features_text)


def write_positions(position_rows: list[list[str]], output_path: Path) -> None:
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(position_rows)


def write_positions_csv(position_rows: list[list[str]], output_path: Path) -> None:
    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        writer.writerows(position_rows)


def write_summary(summary: dict, output_path: Path) -> None:
    with open(output_path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def parse_stat_file(stat_path: Path) -> dict[str, str]:
    stats = {}
    if not stat_path.exists():
        return stats
    with open(stat_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            stats[key.strip()] = value.strip()
    return stats


def write_filtered_stat(
    stat_path: Path,
    original_stats: dict[str, str],
    matrix,
    summary: dict,
) -> None:
    if matrix.shape[1] == 0:
        umi_per_bin = np.array([0])
        genes_per_bin = np.array([0])
        expressed_genes = 0
    else:
        umi_per_bin = np.asarray(matrix.sum(axis=0)).ravel()
        genes_per_bin = np.asarray(matrix.getnnz(axis=0)).ravel()
        expressed_genes = int(np.count_nonzero(matrix.getnnz(axis=1)))

    stats = dict(original_stats)
    stats.update(
        {
            "Total Genes": f"{expressed_genes: ,.0f}",
            "Estimated Number of square bin": f"{int(matrix.shape[1]): ,.0f}",
            "Mean UMI per square bin": f"{int(umi_per_bin.mean()): ,.0f}",
            "Median UMI per square bin": f"{int(np.median(umi_per_bin)): ,.0f}",
            "Mean Genes per square bin": f"{int(genes_per_bin.mean()): ,.0f}",
            "Median Genes per square bin": f"{int(np.median(genes_per_bin)): ,.0f}",
            "Matrix Nonzero Entries": f"{int(matrix.nnz): ,.0f}",
            "Cavity Filter Applied": "true",
            "Cavity Filter Total Barcodes": f"{summary['total_barcodes']: ,.0f}",
            "Cavity Filter Kept Barcodes": f"{summary['kept_barcodes']: ,.0f}",
            "Cavity Filter Removed Barcodes": f"{summary['removed_barcodes']: ,.0f}",
            "Cavity Filter Removed Fraction": f"{summary['removed_fraction'] * 100:.2f}%",
        }
    )

    preferred_order = [
        "Total Genes",
        "Estimated Number of square bin",
        "Fraction Reads in Square",
        "Mean Reads per square bin",
        "Median Reads per square bin",
        "Mean UMI per square bin",
        "Median UMI per square bin",
        "Mean Genes per square bin",
        "Median Genes per square bin",
        "Matrix Nonzero Entries",
        "Cavity Filter Applied",
        "Cavity Filter Total Barcodes",
        "Cavity Filter Kept Barcodes",
        "Cavity Filter Removed Barcodes",
        "Cavity Filter Removed Fraction",
    ]
    written = set()
    with open(stat_path, "w", encoding="utf-8") as handle:
        for key in preferred_order:
            if key in stats:
                handle.write(f"{key}: {stats[key]}\n")
                written.add(key)
        for key, value in stats.items():
            if key not in written:
                handle.write(f"{key}: {value}\n")


def filter_bin_dir(
    bin_dir: str | Path,
    mask_path: str | Path,
    output_dir: str | Path,
    remove_mask_value: str = "black",
    mask_threshold: int = 127,
    dilate: int = 0,
    barcode_radius: int = 0,
    copy_extra_files: bool = True,
    dry_run: bool = False,
    mask_space: str = "auto",
    hires_image: str | Path | None = None,
) -> dict:
    bin_dir = Path(bin_dir).resolve()
    output_dir = Path(output_dir).resolve()
    mask_path = Path(mask_path).resolve()
    matrix_dir = bin_dir / "filtered_feature_bc_matrix"
    spatial_dir = bin_dir / "spatial"

    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    features_path = matrix_dir / "features.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    positions_path = spatial_dir / "tissue_positions_list.csv"
    for required in [barcodes_path, features_path, matrix_path, positions_path, mask_path]:
        if not required.exists():
            raise FileNotFoundError(f"Required input not found: {required}")

    barcodes = load_barcodes(barcodes_path)
    position_rows = load_positions(positions_path)
    position_map = {row[0]: row for row in position_rows}
    missing_positions = [bc for bc in barcodes if bc not in position_map]
    if missing_positions:
        raise ValueError(f"{len(missing_positions)} barcodes missing from tissue_positions_list.csv")

    remove_mask = load_mask(mask_path, threshold=mask_threshold, remove_mask_value=remove_mask_value, dilate=dilate)
    hires_image_path = Path(hires_image).resolve() if hires_image else spatial_dir / "tissue_hires_image.png"
    remove_mask, mask_resize_info = maybe_resize_mask_to_hires(remove_mask, hires_image_path, position_rows, mask_space)

    keep_flags = np.array(
        [not should_remove(position_map[barcode], remove_mask, radius=barcode_radius) for barcode in barcodes],
        dtype=bool,
    )
    total_barcodes = len(barcodes)
    kept_barcodes = int(keep_flags.sum())
    removed_barcodes = total_barcodes - kept_barcodes
    summary = {
        "input_bin_dir": str(bin_dir),
        "output_bin_dir": str(output_dir),
        "mask": str(mask_path),
        "remove_mask_value": remove_mask_value,
        "mask_threshold": int(mask_threshold),
        "dilate": int(dilate),
        "barcode_radius": int(barcode_radius),
        "mask_space": mask_space,
        "hires_image": str(hires_image_path),
        "dry_run": bool(dry_run),
        "total_barcodes": int(total_barcodes),
        "kept_barcodes": kept_barcodes,
        "removed_barcodes": removed_barcodes,
        "removed_fraction": removed_barcodes / total_barcodes if total_barcodes else 0.0,
    }
    if mask_resize_info is not None:
        summary["mask_resize"] = mask_resize_info

    if dry_run:
        return summary

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Output bin directory exists and is not empty: {output_dir}")

    output_matrix_dir = output_dir / "filtered_feature_bc_matrix"
    output_spatial_dir = output_dir / "spatial"
    output_matrix_dir.mkdir(parents=True, exist_ok=True)
    output_spatial_dir.mkdir(parents=True, exist_ok=True)

    features_text = load_features(features_path)
    matrix = sio.mmread(str(matrix_path)).tocsr()
    filtered_matrix = matrix[:, keep_flags].tocsr()
    filtered_barcodes = [barcode for barcode, keep in zip(barcodes, keep_flags) if keep]
    filtered_positions = [position_map[barcode] for barcode in filtered_barcodes]

    tmp_matrix = output_matrix_dir / "matrix.mtx"
    sio.mmwrite(str(tmp_matrix), filtered_matrix.tocoo())
    with open(tmp_matrix, "rb") as src, gzip.open(output_matrix_dir / "matrix.mtx.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    tmp_matrix.unlink()
    write_barcodes(filtered_barcodes, output_matrix_dir / "barcodes.tsv.gz")
    write_features(features_text, output_matrix_dir / "features.tsv.gz")
    write_positions(filtered_positions, output_spatial_dir / "tissue_positions_list.csv")
    write_positions_csv(filtered_positions, output_spatial_dir / "tissue_positions.csv")

    if copy_extra_files:
        for rel_path in [
            Path("downsample.tsv"),
            Path("spatial") / "scalefactors_json.json",
            Path("spatial") / "tissue_lowres_image.png",
        ]:
            src = bin_dir / rel_path
            dst = output_dir / rel_path
            if src.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
        if hires_image_path.exists():
            dst = output_dir / "spatial" / "tissue_hires_image.png"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(hires_image_path, dst)

    write_filtered_stat(output_dir / "stat.txt", parse_stat_file(bin_dir / "stat.txt"), filtered_matrix, summary)
    write_summary(summary, output_dir / "filter_summary.json")
    return summary


def infer_sample_from_square_bin(square_bin_dir: str | Path) -> str:
    square_bin_dir = Path(square_bin_dir)
    for pattern in ("*_bin20", "*_bin10", "*_bin50", "*_bin100", "*_Raw"):
        matches = sorted(square_bin_dir.glob(pattern))
        if matches:
            name = matches[0].name
            for suffix in ("_bin20", "_bin10", "_bin50", "_bin100", "_Raw"):
                if name.endswith(suffix):
                    return name[: -len(suffix)]
    raise ValueError(f"Could not infer sample from square-bin directory: {square_bin_dir}")


def filter_square_bins(
    square_bin_dir: str | Path,
    mask_path: str | Path,
    output_square_bin_dir: str | Path,
    sample: str | None = None,
    bins: tuple[str, ...] = ("10", "20", "50", "100", "Raw"),
    **kwargs,
) -> list[dict]:
    square_bin_dir = Path(square_bin_dir)
    output_square_bin_dir = Path(output_square_bin_dir)
    sample = sample or infer_sample_from_square_bin(square_bin_dir)
    summaries = []
    for bin_label in bins:
        suffix = f"bin{bin_label}" if str(bin_label).lower() != "raw" else "Raw"
        input_bin = square_bin_dir / f"{sample}_{suffix}"
        output_bin = output_square_bin_dir / f"{sample}_{suffix}"
        if not input_bin.exists():
            summaries.append({"bin": suffix, "status": "skipped", "reason": "input_not_found", "input_bin_dir": str(input_bin)})
            continue
        summary = filter_bin_dir(input_bin, mask_path, output_bin, **kwargs)
        summary["bin"] = suffix
        summary["status"] = "dry_run" if kwargs.get("dry_run") else "completed"
        summaries.append(summary)
    return summaries
