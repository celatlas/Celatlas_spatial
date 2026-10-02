#!/usr/bin/env python3
"""
Assign spatial count_detail records to precomputed StarDist label masks.

This script assigns precomputed StarDist labels to spatial barcodes/counts. It
avoids pandas because the current base environment can have a pandas/numpy ABI
mismatch. Large count_detail files are streamed line by line.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import time
from pathlib import Path

import numpy as np
import tifffile


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Assign Celatlas spatial barcodes/counts to an existing StarDist label mask."
    )
    parser.add_argument("--labels", required=True, help="StarDist label TIFF, e.g. *.labels.tif")
    parser.add_argument("--barcode-positions", required=True, help="FilterBarcodes CSV with x,y,barcode rows.")
    parser.add_argument(
        "--barcode-allowlist",
        action="append",
        default=[],
        help=(
            "Optional barcode allowlist file, e.g. filtered square_bin barcodes.tsv.gz. "
            "Can be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--barcode-mask",
        help=(
            "Optional binary tissue mask used to keep barcode positions by coordinate. "
            "Use this to apply cavityFilter tissue masks to raw barcode positions."
        ),
    )
    parser.add_argument(
        "--barcode-mask-keep-value",
        choices=["white", "black"],
        default="white",
        help="Mask pixel value to keep when --barcode-mask is supplied.",
    )
    parser.add_argument("--barcode-mask-threshold", type=int, default=127, help="Threshold for --barcode-mask.")
    parser.add_argument("--barcode-in-tissue-only", action="store_true", help="If barcode positions include in_tissue, keep only in_tissue=1 rows.")
    parser.add_argument("--count-detail", help="05.count/<sample>_count_detail.txt. If omitted, only assignments are written.")
    parser.add_argument("--gtf", help="Reference GTF used to map gene_id to gene_name in features.tsv.gz.")
    parser.add_argument("--tissue-bbox", help="CSV containing bbox as x0,y0,width,height.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--sample", help="Sample name used in output filenames.")
    parser.add_argument("--he-image", help="Registered H&E image used to write StarDist visualization overlays.")
    parser.add_argument("--visualization-preview-dim", type=int, default=1500, help="Long edge for visualization preview images.")
    parser.add_argument(
        "--preview-only",
        action="store_true",
        help="Only regenerate barcode assignment preview images; do not modify assignment or matrix outputs.",
    )
    parser.add_argument("--expand-pixels", type=float, default=0.0, help="Expand nuclear labels by N label-image pixels before assignment.")
    parser.add_argument("--rounding", choices=["floor", "round"], default="floor", help="How to map scaled coordinates to label pixels.")
    parser.add_argument("--chunksize", type=int, default=1_000_000, help="Progress interval for streamed count_detail rows.")
    parser.add_argument("--min-umi", type=int, default=10, help="Minimum assigned UMI rows per retained cell.")
    parser.add_argument("--min-genes", type=int, default=5, help="Minimum detected genes per retained cell.")
    parser.add_argument("--max-umi", type=int, default=None, help="Optional maximum assigned UMI rows per retained cell.")
    parser.add_argument("--no-label-output", action="store_true", help="Do not write expanded/QC relabeled TIFF masks.")
    return parser.parse_args()


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_barcode_allowlist(paths: list[str]) -> set[str] | None:
    if not paths:
        return None
    allowlist: set[str] = set()
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                barcode = line.strip().split("\t", 1)[0].split(",", 1)[0]
                if barcode:
                    allowlist.add(barcode)
    return allowlist


def infer_sample(labels_path: str) -> str:
    name = Path(labels_path).name
    for suffix in (".labels.tif", ".labels.tiff", ".tif", ".tiff"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def parse_bbox(path: str | None, label_shape: tuple[int, int]) -> tuple[float, float, float, float]:
    if not path:
        return (0.0, 0.0, float(label_shape[1]), float(label_shape[0]))

    with open_text(path) as handle:
        for line in handle:
            values = []
            for token in re.split(r"[\t,\s]+", line.strip()):
                if not token:
                    continue
                try:
                    values.append(float(token))
                except ValueError:
                    continue
            if len(values) >= 4:
                x0, y0, width, height = values[:4]
                if width <= 0 or height <= 0:
                    raise ValueError(f"Invalid bbox width/height in {path}: {values[:4]}")
                return x0, y0, width, height

    raise ValueError(f"Failed to parse bbox values from {path}")


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def split_row(line: str) -> list[str]:
    line = line.rstrip("\n\r")
    if "\t" in line and line.count("\t") >= line.count(","):
        return [part.strip() for part in line.split("\t")]
    return [part.strip() for part in line.split(",")]


def read_barcode_positions(
    path: str,
    in_tissue_only: bool = False,
    allowlist: set[str] | None = None,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    barcodes: list[str] = []
    xs: list[float] = []
    ys: list[float] = []
    header: list[str] | None = None

    with open_text(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = split_row(line)
            lowered = [part.lower() for part in parts]
            if line_no == 1 and "barcode" in lowered:
                header = lowered
                continue

            if header:
                try:
                    barcode = parts[header.index("barcode")]
                    if in_tissue_only and "in_tissue" in header:
                        in_tissue_value = parts[header.index("in_tissue")]
                        if in_tissue_value not in {"1", "True", "true", "TRUE"}:
                            continue
                    if {"x", "y"}.issubset(header):
                        x = parts[header.index("x")]
                        y = parts[header.index("y")]
                    elif {"pxl_col_in_fullres", "pxl_row_in_fullres"}.issubset(header):
                        x = parts[header.index("pxl_col_in_fullres")]
                        y = parts[header.index("pxl_row_in_fullres")]
                    else:
                        raise ValueError
                except (ValueError, IndexError) as exc:
                    raise ValueError(f"Unsupported barcode header in {path}: {header}") from exc
            elif len(parts) >= 3 and is_number(parts[0]) and is_number(parts[1]):
                x, y, barcode = parts[0], parts[1], parts[2]
            elif len(parts) >= 6:
                barcode, y, x = parts[0], parts[4], parts[5]
            else:
                raise ValueError(f"Cannot parse barcode row {line_no} in {path}: {parts}")

            if allowlist is not None and str(barcode) not in allowlist:
                continue

            try:
                xs.append(float(x))
                ys.append(float(y))
            except ValueError:
                xs.append(np.nan)
                ys.append(np.nan)
            barcodes.append(str(barcode))

    return barcodes, np.asarray(xs, dtype=np.float64), np.asarray(ys, dtype=np.float64)


def filter_barcode_positions_by_mask(
    barcodes: list[str],
    raw_x: np.ndarray,
    raw_y: np.ndarray,
    mask_path: str,
    bbox: tuple[float, float, float, float],
    keep_value: str = "white",
    threshold: int = 127,
) -> tuple[list[str], np.ndarray, np.ndarray, dict]:
    from PIL import Image

    Image.MAX_IMAGE_PIXELS = None
    mask = np.asarray(Image.open(mask_path).convert("L"), dtype=np.uint8)
    mask_h, mask_w = mask.shape
    x0, y0, width, height = bbox
    mask_x = np.floor((raw_x - x0) * (mask_w / width)).astype(np.int64)
    mask_y = np.floor((raw_y - y0) * (mask_h / height)).astype(np.int64)
    in_bounds = (mask_x >= 0) & (mask_x < mask_w) & (mask_y >= 0) & (mask_y < mask_h)

    foreground = mask > int(threshold)
    if keep_value == "black":
        foreground = ~foreground

    keep = np.zeros(raw_x.shape[0], dtype=bool)
    keep[in_bounds] = foreground[mask_y[in_bounds], mask_x[in_bounds]]
    filtered_indices = np.nonzero(keep)[0]
    summary = {
        "mask": os.path.abspath(mask_path),
        "mask_shape": [int(mask_h), int(mask_w)],
        "keep_value": keep_value,
        "threshold": int(threshold),
        "input_barcodes": int(len(barcodes)),
        "in_mask_bounds_barcodes": int(in_bounds.sum()),
        "kept_barcodes": int(keep.sum()),
        "removed_barcodes": int(len(barcodes) - keep.sum()),
    }
    return (
        [barcodes[idx] for idx in filtered_indices],
        raw_x[filtered_indices],
        raw_y[filtered_indices],
        summary,
    )


def scale_coords(
    raw_x: np.ndarray,
    raw_y: np.ndarray,
    label_shape: tuple[int, int],
    bbox: tuple[float, float, float, float],
    rounding: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x0, y0, width, height = bbox
    label_h, label_w = label_shape
    scaled_x = (raw_x - x0) * (label_w / width)
    scaled_y = (raw_y - y0) * (label_h / height)

    if rounding == "round":
        label_x = np.rint(scaled_x)
        label_y = np.rint(scaled_y)
    else:
        label_x = np.floor(scaled_x)
        label_y = np.floor(scaled_y)

    invalid = ~np.isfinite(label_x) | ~np.isfinite(label_y)
    label_x[invalid] = -1
    label_y[invalid] = -1
    label_x = label_x.astype(np.int64)
    label_y = label_y.astype(np.int64)
    in_bounds = (label_x >= 0) & (label_x < label_w) & (label_y >= 0) & (label_y < label_h)
    return label_x, label_y, in_bounds


def assign_barcodes(labels: np.ndarray, label_x: np.ndarray, label_y: np.ndarray, in_bounds: np.ndarray) -> np.ndarray:
    cell_ids = np.zeros(label_x.shape[0], dtype=np.uint32)
    cell_ids[in_bounds] = labels[label_y[in_bounds], label_x[in_bounds]].astype(np.uint32)
    return cell_ids


def write_boundary_overlay(
    he_image: str,
    qc_labels: np.ndarray,
    output_path: Path,
    preview_path: Path,
    preview_dim: int,
    id_lookup: np.ndarray | None = None,
) -> None:
    from PIL import Image
    from skimage.segmentation import find_boundaries

    Image.MAX_IMAGE_PIXELS = None
    try:
        resample = Image.Resampling.BILINEAR
    except AttributeError:
        resample = Image.BILINEAR

    # A full-resolution 27k x 32k RGB overlay is another ~2.7 GB allocation.
    # For slide-sized inputs, render only a sampled preview from the TIFF
    # memmap and keep the label mask sampled as well.
    if max(qc_labels.shape) > 8000:
        from celatlas_spatial.tools.stardist_tiled_inference import ImageRegionSource

        with ImageRegionSource(Path(he_image)) as source:
            image = source.preview_rgb(max(1, int(preview_dim)))
        sampled = sampled_labels(qc_labels, max(image.size))
        if id_lookup is not None:
            sampled = id_lookup[sampled]
        overlay = np.asarray(image, dtype=np.uint8).copy()
        overlay[find_boundaries(sampled, mode="outer")] = np.array([255, 32, 32], dtype=np.uint8)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        Image.fromarray(overlay).save(output_path, quality=92)
        Image.fromarray(overlay).save(preview_path, quality=92)
        return

    image = Image.open(he_image).convert("RGB").resize((qc_labels.shape[1], qc_labels.shape[0]), resample)
    overlay = np.asarray(image, dtype=np.uint8).copy()
    render_labels = id_lookup[qc_labels] if id_lookup is not None else qc_labels
    boundaries = find_boundaries(render_labels, mode="outer")
    overlay[boundaries] = np.array([255, 32, 32], dtype=np.uint8)
    overlay_image = Image.fromarray(overlay)
    overlay_image.save(output_path, quality=92)
    if preview_dim > 0 and max(overlay_image.size) > preview_dim:
        scale = float(preview_dim) / float(max(overlay_image.size))
        preview_size = (
            max(1, int(round(overlay_image.size[0] * scale))),
            max(1, int(round(overlay_image.size[1] * scale))),
        )
        overlay_image = overlay_image.resize(preview_size, resample)
    overlay_image.save(preview_path, quality=92)


def write_cell_umi_heatmap(
    labels: np.ndarray,
    cell_umi: np.ndarray,
    output_path: Path,
    max_dim: int,
) -> None:
    from PIL import Image

    try:
        import matplotlib.cm as cm
    except Exception:
        cm = None

    try:
        resample = Image.Resampling.BILINEAR
    except AttributeError:
        resample = Image.BILINEAR

    render_labels = sampled_labels(labels, max_dim) if max(labels.shape) > 8000 else labels
    umi_lookup = np.zeros(labels.max() + 1, dtype=np.float32)
    umi_lookup[: min(len(umi_lookup), len(cell_umi))] = cell_umi[: min(len(umi_lookup), len(cell_umi))]
    heat = np.log1p(umi_lookup[render_labels])
    if float(heat.max()) > 0:
        heat = heat / float(heat.max())

    if cm is not None:
        rgba = cm.magma(heat)
        rgba[..., 3] = (render_labels > 0) * 0.88
        rgb = (rgba[..., :3] * 255).astype(np.uint8)
        rgb[render_labels == 0] = 245
        image = Image.fromarray(rgb)
    else:
        image = Image.fromarray((heat * 255).astype(np.uint8)).convert("RGB")

    if max_dim > 0:
        image.thumbnail((max_dim, max_dim), resample)
    image.save(output_path)


def write_high_umi_crops(
    overlay_path: Path,
    sample: str,
    visualization_dir: Path,
    retained: np.ndarray,
    centroid_x: np.ndarray,
    centroid_y: np.ndarray,
    cell_umi: np.ndarray,
    cell_genes: np.ndarray,
    label_shape: tuple[int, int] | None = None,
    crop_size: int = 1000,
    crop_count: int = 4,
) -> None:
    from PIL import Image, ImageDraw

    overlay = Image.open(overlay_path).convert("RGB")
    if label_shape is not None and overlay.size != (label_shape[1], label_shape[0]):
        log("Skipping high-UMI crops because the slide-sized QC overlay is preview-only")
        return
    crop_dir = visualization_dir / "crops"
    crop_dir.mkdir(parents=True, exist_ok=True)

    selected: list[int] = []
    min_dist2 = (crop_size * 0.65) ** 2
    for raw_id in sorted(retained.tolist(), key=lambda idx: int(cell_umi[idx]), reverse=True):
        cx = float(centroid_x[raw_id])
        cy = float(centroid_y[raw_id])
        if not np.isfinite(cx) or not np.isfinite(cy):
            continue
        if all((cx - float(centroid_x[other])) ** 2 + (cy - float(centroid_y[other])) ** 2 > min_dist2 for other in selected):
            selected.append(int(raw_id))
        if len(selected) >= crop_count:
            break

    for idx, raw_id in enumerate(selected, start=1):
        cx = int(round(float(centroid_x[raw_id])))
        cy = int(round(float(centroid_y[raw_id])))
        half = crop_size // 2
        left = max(0, min(overlay.width - crop_size, cx - half))
        top = max(0, min(overlay.height - crop_size, cy - half))
        crop = overlay.crop((left, top, left + crop_size, top + crop_size))
        draw = ImageDraw.Draw(crop, "RGBA")
        draw.rectangle((0, 0, crop.width, 44), fill=(255, 255, 255, 210))
        draw.text(
            (12, 12),
            f"raw={raw_id} UMI={int(cell_umi[raw_id])} genes={int(cell_genes[raw_id])}",
            fill=(0, 0, 0, 255),
        )
        crop.save(crop_dir / f"{sample}.high_umi_crop_{idx:02d}.jpg", quality=94)


def write_barcode_assignment_preview(
    he_image: str,
    raw_x: np.ndarray,
    raw_y: np.ndarray,
    cell_ids: np.ndarray,
    bbox: tuple[float, float, float, float],
    output_path: Path,
    preview_dim: int,
    full_chip_output_path: Path | None = None,
) -> None:
    from PIL import Image, ImageDraw

    Image.MAX_IMAGE_PIXELS = None
    try:
        resample = Image.Resampling.BILINEAR
    except AttributeError:
        resample = Image.BILINEAR

    # PIL crop/resize can decode the complete TIFF. Use the same region-aware
    # reader as tiled inference for very large registered images.
    from celatlas_spatial.tools.stardist_tiled_inference import ImageRegionSource

    with ImageRegionSource(Path(he_image)) as source_region:
        source_width, source_height = source_region.size
        x0, y0, width, height = bbox
        source_x = (raw_x - x0) * (source_width / width)
        source_y = (raw_y - y0) * (source_height / height)
        if max(source_width, source_height) > 8000:
            image = source_region.preview_rgb(max(1, int(preview_dim) or 3000))
            scale_x = image.size[0] / source_width
            scale_y = image.size[1] / source_height
            draw = ImageDraw.Draw(image, "RGBA")
            visible = (
                np.isfinite(source_x)
                & np.isfinite(source_y)
                & (source_x >= 0)
                & (source_x < source_width)
                & (source_y >= 0)
                & (source_y < source_height)
            )
            for x, y in zip(source_x[visible & (cell_ids == 0)], source_y[visible & (cell_ids == 0)]):
                draw.point((int(x * scale_x), int(y * scale_y)), fill=(255, 145, 0, 105))
            for x, y in zip(source_x[visible & (cell_ids > 0)], source_y[visible & (cell_ids > 0)]):
                px, py = int(x * scale_x), int(y * scale_y)
                draw.ellipse((px - 1, py - 1, px + 1, py + 1), fill=(0, 220, 255, 155))
            output_path.parent.mkdir(parents=True, exist_ok=True)
            image.save(output_path, quality=92)
            if full_chip_output_path is not None:
                full_chip_output_path.parent.mkdir(parents=True, exist_ok=True)
                image.save(full_chip_output_path, quality=92)
            return

    source = Image.open(he_image).convert("RGB")
    source_x = (raw_x - x0) * (source.size[0] / width)
    source_y = (raw_y - y0) * (source.size[1] / height)
    in_bounds = (
        np.isfinite(source_x)
        & np.isfinite(source_y)
        & (source_x >= 0)
        & (source_x < source.size[0])
        & (source_y >= 0)
        & (source_y < source.size[1])
    )
    assigned = in_bounds & (cell_ids > 0)

    def expanded_interval(values: np.ndarray, limit: int) -> tuple[int, int]:
        if values.size == 0:
            return 0, limit
        margin = max(2.0, limit * 0.02)
        low = float(np.nanmin(values)) - margin
        high = float(np.nanmax(values)) + margin
        minimum_span = max(4.0, limit * 0.10)
        if high - low < minimum_span:
            center = (low + high) / 2.0
            low = center - minimum_span / 2.0
            high = center + minimum_span / 2.0
        if low < 0:
            high -= low
            low = 0.0
        if high > limit:
            low -= high - limit
            high = float(limit)
        return max(0, int(np.floor(low))), min(limit, int(np.ceil(high)))

    def render(crop_box: tuple[int, int, int, int], include_unassigned: bool) -> Image.Image:
        image = source.crop(crop_box)
        crop_width = crop_box[2] - crop_box[0]
        crop_height = crop_box[3] - crop_box[1]
        if preview_dim > 0:
            image.thumbnail((preview_dim, preview_dim), resample)
        x_scale = image.size[0] / crop_width
        y_scale = image.size[1] / crop_height
        xs = np.floor((source_x - crop_box[0]) * x_scale).astype(np.int64)
        ys = np.floor((source_y - crop_box[1]) * y_scale).astype(np.int64)
        visible = (
            in_bounds
            & (xs >= 0)
            & (xs < image.size[0])
            & (ys >= 0)
            & (ys < image.size[1])
        )
        draw = ImageDraw.Draw(image, "RGBA")
        if include_unassigned:
            unassigned = visible & (cell_ids == 0)
            # Unassigned barcodes are dense; one-pixel points preserve the HE background.
            for x, y in zip(xs[unassigned], ys[unassigned]):
                draw.point((int(x), int(y)), fill=(255, 145, 0, 105))
        visible_assigned = visible & (cell_ids > 0)
        for x, y in zip(xs[visible_assigned], ys[visible_assigned]):
            draw.ellipse(
                (int(x) - 1, int(y) - 1, int(x) + 1, int(y) + 1),
                fill=(0, 220, 255, 155),
            )
        return image

    if assigned.any():
        left, right = expanded_interval(source_x[assigned], source.size[0])
        top, bottom = expanded_interval(source_y[assigned], source.size[1])
        tissue_crop = (left, top, right, bottom)
    else:
        tissue_crop = (0, 0, source.size[0], source.size[1])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    render(tissue_crop, include_unassigned=False).save(output_path, quality=92)
    if full_chip_output_path is not None:
        full_chip_output_path.parent.mkdir(parents=True, exist_ok=True)
        render((0, 0, source.size[0], source.size[1]), include_unassigned=True).save(
            full_chip_output_path,
            quality=92,
        )


def write_assignments(
    path: Path,
    barcodes: list[str],
    raw_x: np.ndarray,
    raw_y: np.ndarray,
    label_x: np.ndarray,
    label_y: np.ndarray,
    in_bounds: np.ndarray,
    cell_ids: np.ndarray,
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("barcode\traw_x\traw_y\tlabel_x\tlabel_y\tin_bounds\tcell_id\n")
        for idx, barcode in enumerate(barcodes):
            handle.write(
                f"{barcode}\t{raw_x[idx]:.6f}\t{raw_y[idx]:.6f}\t"
                f"{int(label_x[idx])}\t{int(label_y[idx])}\t{int(in_bounds[idx])}\t{int(cell_ids[idx])}\n"
            )


def write_per_cell_barcodes(path: Path, cell_ids: np.ndarray, max_label: int) -> np.ndarray:
    assigned = cell_ids[cell_ids > 0].astype(np.int64)
    counts = np.bincount(assigned, minlength=max_label + 1).astype(np.int64)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("raw_cell_id\tassigned_barcodes\n")
        for cell_id in np.nonzero(counts)[0]:
            if cell_id == 0:
                continue
            handle.write(f"{int(cell_id)}\t{int(counts[cell_id])}\n")
    return counts


def assignment_summary(cell_ids: np.ndarray, in_bounds: np.ndarray, max_label: int) -> dict:
    total = int(cell_ids.shape[0])
    assigned = cell_ids[cell_ids > 0].astype(np.int64)
    per_cell = np.bincount(assigned, minlength=max_label + 1) if assigned.size else np.zeros(max_label + 1)
    nonzero = per_cell[per_cell > 0]
    if nonzero.size:
        quantiles = np.quantile(nonzero, [0.1, 0.25, 0.5, 0.75, 0.9])
        quantile_dict = {
            "0.10": float(quantiles[0]),
            "0.25": float(quantiles[1]),
            "0.50": float(quantiles[2]),
            "0.75": float(quantiles[3]),
            "0.90": float(quantiles[4]),
        }
    else:
        quantile_dict = {}
    return {
        "total_labels": int(max_label),
        "total_barcodes": total,
        "in_bounds_barcodes": int(in_bounds.sum()),
        "assigned_barcodes": int(assigned.size),
        "assignment_rate": round(float(assigned.size / total), 6) if total else 0.0,
        "assigned_cells": int(nonzero.size),
        "assigned_cell_fraction": round(float(nonzero.size / max_label), 6) if max_label else 0.0,
        "barcodes_per_assigned_cell": quantile_dict,
    }


def write_json(path: Path, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def parse_count_header(first_line: str) -> tuple[dict[str, int], list[str] | None]:
    parts = first_line.rstrip("\n\r").split("\t")
    lowered = [part.strip().lower() for part in parts]
    synonyms = {
        "barcode": {"barcode"},
        "gene": {"geneid", "gene_id", "gene"},
        "umi": {"umi"},
        "count": {"count", "read_count"},
    }
    indices: dict[str, int] = {}
    for name, aliases in synonyms.items():
        for idx, col in enumerate(lowered):
            if col in aliases:
                indices[name] = idx
                break

    if {"barcode", "gene", "umi"}.issubset(indices):
        return indices, None

    return {"barcode": 0, "gene": 1, "umi": 2, "count": 3}, parts


def encode_pair(gene_idx: int, cell_id: int) -> int:
    return (int(gene_idx) << 32) | int(cell_id)


def decode_pair(key: int) -> tuple[int, int]:
    return key >> 32, key & 0xFFFFFFFF


def parse_gene_id_name_map(gtf_path: str | None) -> dict[str, str]:
    if not gtf_path:
        return {}
    path = Path(gtf_path)
    if not path.exists():
        raise FileNotFoundError(f"GTF not found: {gtf_path}")

    opener = gzip.open if str(path).endswith(".gz") else open
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')
    gene_name_pattern = re.compile(r'gene_name "([^"]+)"')
    gene_map: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gene_id_match = gene_id_pattern.search(parts[8])
            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)
            gene_name_match = gene_name_pattern.search(parts[8])
            gene_name = gene_name_match.group(1) if gene_name_match else gene_id
            gene_map[gene_id] = gene_name
    return gene_map


def write_10x_matrix(
    matrix,
    genes: list[str],
    barcodes: list[str],
    outdir: Path,
    gene_id_name_map: dict[str, str] | None = None,
) -> None:
    import scipy.io

    gene_id_name_map = gene_id_name_map or {}
    outdir.mkdir(parents=True, exist_ok=True)
    with gzip.open(outdir / "features.tsv.gz", "wt", encoding="utf-8") as handle:
        for gene in genes:
            gene_name = gene_id_name_map.get(gene, gene)
            handle.write(f"{gene}\t{gene_name}\tGene Expression\n")

    with gzip.open(outdir / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
        for barcode in barcodes:
            handle.write(f"{barcode}\n")

    tmp_mtx = outdir / "matrix.mtx"
    scipy.io.mmwrite(str(tmp_mtx), matrix.tocoo())
    with open(tmp_mtx, "rb") as source, gzip.open(outdir / "matrix.mtx.gz", "wb") as target:
        shutil.copyfileobj(source, target)
    tmp_mtx.unlink()


def compute_label_geometry(
    labels: np.ndarray,
    max_label: int,
    bbox: tuple[float, float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    label_h, label_w = labels.shape
    area = np.zeros(max_label + 1, dtype=np.int64)
    sum_x = np.zeros(max_label + 1, dtype=np.float64)
    sum_y = np.zeros(max_label + 1, dtype=np.float64)
    x_coords = np.arange(label_w, dtype=np.float64)

    for row_idx in range(label_h):
        row = labels[row_idx]
        mask = row > 0
        if not mask.any():
            continue
        ids = row[mask].astype(np.int64, copy=False)
        area += np.bincount(ids, minlength=max_label + 1)
        sum_y += np.bincount(ids, weights=np.full(ids.shape, row_idx, dtype=np.float64), minlength=max_label + 1)
        sum_x += np.bincount(ids, weights=x_coords[mask], minlength=max_label + 1)

    valid = area > 0
    centroid_x = np.zeros(max_label + 1, dtype=np.float64)
    centroid_y = np.zeros(max_label + 1, dtype=np.float64)
    centroid_x[valid] = sum_x[valid] / area[valid]
    centroid_y[valid] = sum_y[valid] / area[valid]

    x0, y0, width, height = bbox
    fullres_x = x0 + (centroid_x + 0.5) * (width / label_w)
    fullres_y = y0 + (centroid_y + 0.5) * (height / label_h)
    return area, centroid_x, centroid_y, fullres_x, fullres_y


def open_label_memmap(path: str) -> np.ndarray:
    """Open an uncompressed label TIFF without copying its full pixel plane."""
    try:
        labels = tifffile.memmap(path, mode="r")
        if labels.ndim != 2:
            raise ValueError(f"Expected a 2D label mask, got shape {labels.shape}")
        return labels
    except (OSError, ValueError, tifffile.TiffFileError):
        labels = tifffile.imread(path)
        if labels.ndim != 2:
            raise ValueError(f"Expected a 2D label mask, got shape {labels.shape}")
        return labels.astype(np.uint32, copy=False)


def expand_labels_chunked(
    labels: np.ndarray,
    distance: float,
    output_path: Path,
    chunk_rows: int = 1024,
) -> np.memmap:
    """Expand a slide-sized label mask in slabs, keeping memory bounded.

    Each slab includes a halo wider than the requested expansion distance, so
    pixels written to the slab interior see the same neighbourhood as a
    whole-image expansion.  The output is an uncompressed uint32 BigTIFF
    memmap and can therefore be consumed by the streaming assignment code.
    """
    from skimage.segmentation import expand_labels

    if distance <= 0:
        raise ValueError("distance must be positive")
    radius = max(1, int(np.ceil(float(distance))))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    expanded = tifffile.memmap(
        output_path,
        shape=labels.shape,
        dtype=np.uint32,
        bigtiff=True,
        photometric="minisblack",
    )
    height = labels.shape[0]
    for y0 in range(0, height, max(1, int(chunk_rows))):
        y1 = min(height, y0 + max(1, int(chunk_rows)))
        read_y0 = max(0, y0 - radius)
        read_y1 = min(height, y1 + radius)
        slab = np.asarray(labels[read_y0:read_y1], dtype=np.uint32)
        slab_expanded = expand_labels(slab, distance=float(distance)).astype(np.uint32, copy=False)
        expanded[y0:y1] = slab_expanded[y0 - read_y0 : y1 - read_y0]
        del slab, slab_expanded
    expanded.flush()
    return expanded


def sampled_labels(labels: np.ndarray, max_dim: int) -> np.ndarray:
    """Return a small label preview without loading the full label canvas."""
    height, width = labels.shape
    scale = min(1.0, float(max_dim) / max(width, height)) if max_dim > 0 else 1.0
    preview_height = max(1, round(height * scale))
    preview_width = max(1, round(width * scale))
    y_index = np.minimum(
        (np.arange(preview_height, dtype=np.float64) / scale).astype(np.int64), height - 1
    )
    x_index = np.minimum(
        (np.arange(preview_width, dtype=np.float64) / scale).astype(np.int64), width - 1
    )
    return np.asarray(labels[np.ix_(y_index, x_index)], dtype=np.uint32)


def remap_labels_chunked(
    labels: np.ndarray,
    id_lookup: np.ndarray,
    output_path: Path,
    chunk_rows: int = 1024,
) -> np.memmap:
    """Write a compacted uint32 label TIFF without a full-size RAM copy."""
    output = tifffile.memmap(
        output_path,
        shape=labels.shape,
        dtype=np.uint32,
        bigtiff=True,
        photometric="minisblack",
    )
    for y0 in range(0, labels.shape[0], max(1, int(chunk_rows))):
        y1 = min(labels.shape[0], y0 + max(1, int(chunk_rows)))
        output[y0:y1] = id_lookup[np.asarray(labels[y0:y1], dtype=np.uint32)]
    output.flush()
    return output


def write_count_outputs(
    outdir: Path,
    sample: str,
    labels: np.ndarray,
    bbox: tuple[float, float, float, float],
    pair_counts: dict[int, int],
    gene_to_idx: dict[str, int],
    cell_umi: np.ndarray,
    cell_reads: np.ndarray,
    assigned_barcode_counts: np.ndarray,
    min_umi: int,
    min_genes: int,
    max_umi: int | None,
    write_labels: bool,
    he_image: str | None,
    visualization_preview_dim: int,
    gene_id_name_map: dict[str, str] | None = None,
) -> dict:
    import scipy.sparse

    max_label = int(labels.max())
    idx_to_gene = [None] * len(gene_to_idx)
    for gene, idx in gene_to_idx.items():
        idx_to_gene[idx] = gene
    gene_order = sorted(range(len(idx_to_gene)), key=lambda idx: idx_to_gene[idx])
    genes_sorted = [idx_to_gene[idx] for idx in gene_order]
    row_remap = np.zeros(len(idx_to_gene), dtype=np.int32)
    for new_idx, old_idx in enumerate(gene_order):
        row_remap[old_idx] = new_idx

    rows = np.empty(len(pair_counts), dtype=np.int32)
    cols = np.empty(len(pair_counts), dtype=np.int32)
    data = np.empty(len(pair_counts), dtype=np.int32)
    cell_genes = np.zeros(max_label + 1, dtype=np.int64)
    for out_idx, (key, value) in enumerate(pair_counts.items()):
        gene_idx, cell_id = decode_pair(key)
        rows[out_idx] = row_remap[gene_idx]
        cols[out_idx] = cell_id
        data[out_idx] = value
        cell_genes[cell_id] += 1

    raw_matrix = scipy.sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(len(genes_sorted), max_label + 1),
        dtype=np.int32,
    ).tocsc()

    keep = (cell_umi >= min_umi) & (cell_genes >= min_genes)
    keep[0] = False
    if max_umi is not None:
        keep &= cell_umi <= max_umi
    retained = np.nonzero(keep)[0].astype(np.int64)
    if retained.size == 0:
        raise ValueError("No cells retained after QC. Lower --min-umi/--min-genes or inspect barcode alignment.")

    barcodes = [f"cell_{idx}" for idx in range(1, retained.size + 1)]
    final_matrix = raw_matrix[:, retained].tocoo()
    matrix_dir = outdir / "cell_matrix"
    write_10x_matrix(final_matrix, genes_sorted, barcodes, matrix_dir, gene_id_name_map=gene_id_name_map)

    new_id_lookup = np.zeros(max_label + 1, dtype=np.int64)
    new_id_lookup[retained] = np.arange(1, retained.size + 1, dtype=np.int64)

    with open(outdir / f"{sample}.cell_id_map.tsv", "w", encoding="utf-8") as handle:
        handle.write("raw_cell_id\tnew_cell_id\tBarcode\n")
        for raw_id in retained:
            new_id = int(new_id_lookup[raw_id])
            handle.write(f"{int(raw_id)}\t{new_id}\tcell_{new_id}\n")

    with open(outdir / f"{sample}.cell_stats.tsv", "w", encoding="utf-8") as handle:
        handle.write("raw_cell_id\tnew_cell_id\tBarcode\tretained\tumi\tgenes\treads\tassigned_barcodes\n")
        for raw_id in range(1, max_label + 1):
            new_id = int(new_id_lookup[raw_id])
            barcode = f"cell_{new_id}" if new_id else ""
            retained_flag = 1 if new_id else 0
            handle.write(
                f"{raw_id}\t{new_id}\t{barcode}\t{retained_flag}\t"
                f"{int(cell_umi[raw_id])}\t{int(cell_genes[raw_id])}\t"
                f"{int(cell_reads[raw_id])}\t{int(assigned_barcode_counts[raw_id])}\n"
            )

    log("Computing label geometry for retained-cell metadata")
    area, centroid_x, centroid_y, fullres_x, fullres_y = compute_label_geometry(labels, max_label, bbox)
    with open(outdir / f"{sample}.cell_metadata.tsv", "w", encoding="utf-8") as handle:
        handle.write(
            "Barcode\traw_cell_id\tnew_cell_id\tarea_px\tcentroid_col_px\tcentroid_row_px\t"
            "pxl_col_in_fullres\tpxl_row_in_fullres\tumi\tgenes\treads\tassigned_barcodes\n"
        )
        for raw_id in retained:
            new_id = int(new_id_lookup[raw_id])
            handle.write(
                f"cell_{new_id}\t{int(raw_id)}\t{new_id}\t{int(area[raw_id])}\t"
                f"{centroid_x[raw_id]:.6f}\t{centroid_y[raw_id]:.6f}\t"
                f"{fullres_x[raw_id]:.6f}\t{fullres_y[raw_id]:.6f}\t"
                f"{int(cell_umi[raw_id])}\t{int(cell_genes[raw_id])}\t"
                f"{int(cell_reads[raw_id])}\t{int(assigned_barcode_counts[raw_id])}\n"
            )

    # Build a retained-cell label view for QC visualizations even when the
    # caller deliberately suppresses the extra full-resolution QC label TIFF.
    # Previously visualization generation lived inside ``if write_labels``;
    # --no-label-output therefore also (incorrectly) removed both QC images.
    if write_labels or he_image:
        id_lookup = np.zeros(max_label + 1, dtype=np.uint32)
        id_lookup[retained] = np.arange(1, retained.size + 1, dtype=np.uint32)
        if write_labels:
            qc_labels_path = outdir / f"{sample}.stardist_qc_labels.tif"
            if max(labels.shape) > 8000:
                qc_labels = remap_labels_chunked(labels, id_lookup, qc_labels_path)
            else:
                qc_labels = id_lookup[labels]
                tifffile.imwrite(qc_labels_path, qc_labels.astype(np.uint32))
        else:
            # Remap only the sampled label pixels in the renderer: no extra
            # full-canvas label array or temporary QC TIFF is needed.
            qc_labels = labels
        if he_image:
            visualization_dir = outdir / "visualization"
            visualization_dir.mkdir(parents=True, exist_ok=True)
            write_boundary_overlay(
                he_image=he_image,
                qc_labels=qc_labels,
                output_path=visualization_dir / f"{sample}.qc_boundary_overlay.jpg",
                preview_path=visualization_dir / f"{sample}.qc_boundary_overlay.preview.jpg",
                preview_dim=visualization_preview_dim,
                id_lookup=None if write_labels else id_lookup,
            )
            write_cell_umi_heatmap(
                labels=labels,
                cell_umi=cell_umi,
                output_path=visualization_dir / f"{sample}.cell_umi_heatmap.png",
                max_dim=2000,
            )
            write_high_umi_crops(
                overlay_path=visualization_dir / f"{sample}.qc_boundary_overlay.jpg",
                sample=sample,
                visualization_dir=visualization_dir,
                retained=retained,
                centroid_x=centroid_x,
                centroid_y=centroid_y,
                cell_umi=cell_umi,
                cell_genes=cell_genes,
                label_shape=labels.shape,
            )
        if isinstance(qc_labels, np.memmap):
            del qc_labels

    retained_umi = cell_umi[retained]
    retained_genes = cell_genes[retained]
    return {
        "genes": int(len(genes_sorted)),
        "matrix_nonzero_entries": int(final_matrix.nnz),
        "cells_before_qc": int(max_label),
        "cells_with_umi": int(np.count_nonzero(cell_umi[1:] > 0)),
        "cells_after_qc": int(retained.size),
        "min_umi": int(min_umi),
        "min_genes": int(min_genes),
        "max_umi": max_umi,
        "median_umi_per_retained_cell": float(np.median(retained_umi)),
        "median_genes_per_retained_cell": float(np.median(retained_genes)),
        "matrix_dir": str(matrix_dir),
    }


def aggregate_counts(
    count_detail: str,
    barcode_to_cell: dict[str, int],
    labels: np.ndarray,
    sample: str,
    outdir: Path,
    bbox: tuple[float, float, float, float],
    assigned_barcode_counts: np.ndarray,
    progress_every: int,
    min_umi: int,
    min_genes: int,
    max_umi: int | None,
    write_labels: bool,
    he_image: str | None,
    visualization_preview_dim: int,
    gtf: str | None = None,
) -> dict:
    max_label = int(labels.max())
    gene_to_idx: dict[str, int] = {}
    pair_counts: dict[int, int] = {}
    cell_umi = np.zeros(max_label + 1, dtype=np.int64)
    cell_reads = np.zeros(max_label + 1, dtype=np.int64)
    total_rows = 0
    assigned_rows = 0
    malformed_rows = 0
    duplicate_barcodes = len(barcode_to_cell) - len(set(barcode_to_cell))
    start = time.time()
    gene_id_name_map = parse_gene_id_name_map(gtf)
    if gene_id_name_map:
        log(f"Loaded gene_id->gene_name map: {len(gene_id_name_map):,} genes")

    with open_text(count_detail) as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError(f"Empty count_detail file: {count_detail}")
        indices, first_data = parse_count_header(first_line)
        max_idx = max(indices.values())

        def process_parts(parts: list[str]) -> None:
            nonlocal assigned_rows, malformed_rows
            if len(parts) <= max_idx:
                malformed_rows += 1
                return
            cell_id = barcode_to_cell.get(parts[indices["barcode"]])
            if not cell_id:
                return
            gene = parts[indices["gene"]]
            gene_idx = gene_to_idx.get(gene)
            if gene_idx is None:
                gene_idx = len(gene_to_idx)
                gene_to_idx[gene] = gene_idx
            count_idx = indices.get("count")
            if count_idx is None or count_idx >= len(parts):
                read_count = 1
            else:
                try:
                    read_count = int(float(parts[count_idx]))
                except ValueError:
                    read_count = 1
            assigned_rows += 1
            cell_umi[cell_id] += 1
            cell_reads[cell_id] += read_count
            key = encode_pair(gene_idx, cell_id)
            pair_counts[key] = pair_counts.get(key, 0) + 1

        if first_data is not None:
            total_rows += 1
            process_parts(first_data)

        log(f"Streaming count_detail with progress interval {progress_every:,} rows")
        for line in handle:
            total_rows += 1
            process_parts(line.rstrip("\n\r").split("\t"))
            if progress_every > 0 and total_rows % progress_every == 0:
                elapsed = (time.time() - start) / 60
                log(
                    f"rows={total_rows:,}, assigned_rows={assigned_rows:,}, "
                    f"genes={len(gene_to_idx):,}, cell_gene_pairs={len(pair_counts):,}, elapsed={elapsed:.1f} min"
                )

    if not pair_counts:
        raise ValueError("No count_detail rows matched barcodes assigned to StarDist labels.")

    summary = write_count_outputs(
        outdir=outdir,
        sample=sample,
        labels=labels,
        bbox=bbox,
        pair_counts=pair_counts,
        gene_to_idx=gene_to_idx,
        cell_umi=cell_umi,
        cell_reads=cell_reads,
        assigned_barcode_counts=assigned_barcode_counts,
        min_umi=min_umi,
        min_genes=min_genes,
        max_umi=max_umi,
        write_labels=write_labels,
        he_image=he_image,
        visualization_preview_dim=visualization_preview_dim,
        gene_id_name_map=gene_id_name_map,
    )
    summary.update(
        {
            "count_detail_rows": int(total_rows),
            "assigned_count_detail_rows": int(assigned_rows),
            "malformed_count_detail_rows": int(malformed_rows),
            "duplicate_assigned_barcodes_dropped": int(duplicate_barcodes),
        }
    )
    return summary


def main() -> None:
    args = parse_args()
    if args.preview_only and not args.he_image:
        raise ValueError("--preview-only requires --he-image")
    sample = args.sample or infer_sample(args.labels)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"Opening labels (memory-mapped when possible): {args.labels}")
    labels = open_label_memmap(args.labels)
    labels = labels.astype(np.uint32, copy=False)
    raw_label_max = int(labels.max())
    bbox = parse_bbox(args.tissue_bbox, labels.shape)
    log(f"label_shape={labels.shape}, raw_labels={raw_label_max:,}, bbox={bbox}")

    transient_expanded_path: Path | None = None
    if args.expand_pixels > 0:
        log(f"Expanding labels by {args.expand_pixels} label pixels")
        if args.no_label_output:
            transient_expanded_path = outdir / f".{sample}.{os.getpid()}.stardist_expanded.work.tif"
            expanded_path = transient_expanded_path
        else:
            expanded_path = outdir / f"{sample}.stardist_expanded_labels.tif"
        expanded = expand_labels_chunked(labels, float(args.expand_pixels), expanded_path)
        if isinstance(labels, np.memmap):
            del labels
        labels = expanded

    max_label = int(labels.max())
    barcode_allowlist = read_barcode_allowlist(args.barcode_allowlist)
    if barcode_allowlist is not None:
        log(f"Loaded barcode allowlist: {len(barcode_allowlist):,} barcodes")
    barcodes, raw_x, raw_y = read_barcode_positions(
        args.barcode_positions,
        args.barcode_in_tissue_only,
        allowlist=barcode_allowlist,
    )
    barcode_mask_summary = None
    if args.barcode_mask:
        barcodes, raw_x, raw_y, barcode_mask_summary = filter_barcode_positions_by_mask(
            barcodes=barcodes,
            raw_x=raw_x,
            raw_y=raw_y,
            mask_path=args.barcode_mask,
            bbox=bbox,
            keep_value=args.barcode_mask_keep_value,
            threshold=args.barcode_mask_threshold,
        )
        log(
            "barcode mask: "
            f"kept {barcode_mask_summary['kept_barcodes']:,}/"
            f"{barcode_mask_summary['input_barcodes']:,} barcodes"
        )
    if not barcodes:
        raise ValueError("No barcode positions remain after allowlist/mask filtering.")

    label_x, label_y, in_bounds = scale_coords(raw_x, raw_y, labels.shape, bbox, args.rounding)
    cell_ids = assign_barcodes(labels, label_x, label_y, in_bounds)

    if args.he_image:
        visualization_dir = outdir / "visualization"
        visualization_dir.mkdir(parents=True, exist_ok=True)
        write_barcode_assignment_preview(
            he_image=args.he_image,
            raw_x=raw_x,
            raw_y=raw_y,
            cell_ids=cell_ids,
            bbox=bbox,
            output_path=visualization_dir / f"{sample}.barcode_assignment.preview.jpg",
            preview_dim=int(args.visualization_preview_dim),
            full_chip_output_path=(
                visualization_dir / f"{sample}.barcode_assignment.full_chip_qc.jpg"
            ),
        )
    if args.preview_only:
        log(f"Barcode assignment preview regenerated: {visualization_dir}")
        if transient_expanded_path is not None:
            transient_expanded_path.unlink(missing_ok=True)
        return

    assignments_path = outdir / f"{sample}.stardist_barcode_assignments.tsv"
    write_assignments(assignments_path, barcodes, raw_x, raw_y, label_x, label_y, in_bounds, cell_ids)
    assigned_barcode_counts = write_per_cell_barcodes(
        outdir / f"{sample}.barcodes_per_cell.tsv",
        cell_ids,
        max_label,
    )

    summary = {
        "sample": sample,
        "labels": os.path.abspath(args.labels),
        "barcode_positions": os.path.abspath(args.barcode_positions),
        "barcode_allowlist": [os.path.abspath(path) for path in args.barcode_allowlist],
        "barcode_allowlist_size": len(barcode_allowlist) if barcode_allowlist is not None else None,
        "barcode_mask": barcode_mask_summary,
        "barcode_in_tissue_only": bool(args.barcode_in_tissue_only),
        "count_detail": os.path.abspath(args.count_detail) if args.count_detail else None,
        "tissue_bbox": os.path.abspath(args.tissue_bbox) if args.tissue_bbox else None,
        "label_shape": list(labels.shape),
        "raw_labels_before_expansion": raw_label_max,
        "labels_after_expansion": max_label,
        "expand_pixels": float(args.expand_pixels),
        "rounding": args.rounding,
        "bbox": list(bbox),
        "assignments_file": str(assignments_path),
        "assignment": assignment_summary(cell_ids, in_bounds, max_label),
    }
    write_json(outdir / f"{sample}.stardist_assignment_summary.json", summary)
    log(
        "assignment: "
        f"{summary['assignment']['assigned_barcodes']:,}/{summary['assignment']['total_barcodes']:,} "
        f"barcodes ({summary['assignment']['assignment_rate'] * 100:.2f}%), "
        f"{summary['assignment']['assigned_cells']:,} cells hit"
    )

    if args.count_detail:
        barcode_to_cell: dict[str, int] = {}
        duplicate_barcodes = 0
        for barcode, cell_id in zip(barcodes, cell_ids):
            if cell_id <= 0:
                continue
            if barcode in barcode_to_cell:
                duplicate_barcodes += 1
                continue
            barcode_to_cell[barcode] = int(cell_id)
        if duplicate_barcodes:
            log(f"WARNING: dropped {duplicate_barcodes:,} duplicate assigned barcodes")

        count_summary = aggregate_counts(
            count_detail=args.count_detail,
            barcode_to_cell=barcode_to_cell,
            labels=labels,
            sample=sample,
            outdir=outdir,
            bbox=bbox,
            assigned_barcode_counts=assigned_barcode_counts,
            progress_every=args.chunksize,
            min_umi=args.min_umi,
            min_genes=args.min_genes,
            max_umi=args.max_umi,
            write_labels=not args.no_label_output,
            he_image=args.he_image,
            visualization_preview_dim=int(args.visualization_preview_dim),
            gtf=args.gtf,
        )
        summary["count_aggregation"] = count_summary
        write_json(outdir / f"{sample}.stardist_assignment_summary.json", summary)
        log(
            "count aggregation: "
            f"{count_summary['cells_after_qc']:,} retained cells, "
            f"{count_summary['genes']:,} genes, matrix={count_summary['matrix_dir']}"
        )
        if args.he_image:
            visualization_dir = outdir / "visualization"
            visualization_dir.mkdir(parents=True, exist_ok=True)
            index_path = visualization_dir / "index.html"
            crop_sections = [
                f"<section><h3>{path.name}</h3><img src=\"crops/{path.name}\"></section>"
                for path in sorted((visualization_dir / "crops").glob(f"{sample}.high_umi_crop_*.jpg"))
            ]
            index_path.write_text(
                "\n".join(
                    [
                        "<!doctype html>",
                        "<html><head><meta charset=\"utf-8\">",
                        f"<title>{sample} StarDist Visualization</title>",
                        (
                            "<style>body{font-family:sans-serif;margin:24px;background:#f7f5ef;color:#1e1e1e} "
                            "img{max-width:100%;border:1px solid #ddd} "
                            ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px}</style>"
                        ),
                        "</head><body>",
                        f"<h1>{sample} StarDist Assignment Visualization</h1>",
                        (
                            f"<p>QC cells: {count_summary['cells_after_qc']:,}; "
                            f"genes: {count_summary['genes']:,}; "
                            f"assigned barcodes: {summary['assignment']['assigned_barcodes']:,}/"
                            f"{summary['assignment']['total_barcodes']:,} "
                            f"({summary['assignment']['assignment_rate']:.2%}).</p>"
                        ),
                        "<div class=\"grid\">",
                        (
                            f"<section><h2>QC Boundary Overlay</h2>"
                            f"<a href=\"{sample}.qc_boundary_overlay.jpg\">"
                            f"<img src=\"{sample}.qc_boundary_overlay.preview.jpg\"></a></section>"
                        ),
                        (
                            f"<section><h2>Barcode Assignment: Tissue View</h2>"
                            f"<img src=\"{sample}.barcode_assignment.preview.jpg\">"
                            f"<p>cyan assigned barcodes.</p></section>"
                        ),
                        (
                            f"<section><h2>Barcode Assignment: Full-chip QC</h2>"
                            f"<img src=\"{sample}.barcode_assignment.full_chip_qc.jpg\">"
                            f"<p>cyan assigned; orange unassigned.</p></section>"
                        ),
                        (
                            f"<section><h2>Cell UMI Heatmap</h2>"
                            f"<img src=\"{sample}.cell_umi_heatmap.png\"></section>"
                        ),
                        "</div>",
                        "<h2>High-UMI Crops</h2>",
                        "<div class=\"grid\">",
                        *crop_sections,
                        "</div></body></html>",
                    ]
                )
            )
            log(f"visualization={index_path}")

    if transient_expanded_path is not None:
        del labels
        transient_expanded_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
