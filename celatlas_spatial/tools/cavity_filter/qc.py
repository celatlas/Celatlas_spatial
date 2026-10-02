"""QC plots for cavity-filtered square-bin matrices."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import io as sio

Image.MAX_IMAGE_PIXELS = None


_CAVITY_OVERLAY_PREFIXES = ("3a_", "3b_", "3c_", "4_", "5_")
_JET_ZERO_RGB = np.array([0, 0, 127], dtype=np.uint8)


def _load_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def _load_positions(path: Path) -> dict[str, tuple[float, float]]:
    positions = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            if row[0] == "barcode":
                continue
            positions[row[0]] = (float(row[4]), float(row[5]))
    return positions


def _image_shape(path: Path) -> tuple[int, int] | None:
    if not path.exists():
        return None
    try:
        with Image.open(path) as image:
            return image.size[1], image.size[0]
    except Exception:
        return None


def _resize_shape(shape: tuple[int, int], max_dim: int) -> tuple[int, int]:
    height, width = shape
    if max_dim <= 0 or max(height, width) <= max_dim:
        return height, width
    scale = float(max_dim) / float(max(height, width))
    return max(1, int(round(height * scale))), max(1, int(round(width * scale)))


def _resize_image(image: Image.Image, max_dim: int, resample: int) -> Image.Image:
    target_height, target_width = _resize_shape((image.height, image.width), max_dim=max_dim)
    if (target_width, target_height) == image.size:
        return image
    return image.resize((target_width, target_height), resample=resample)


def write_cavity_overlay_qc(
    source_overlay_dir: str | Path,
    tissue_mask_path: str | Path,
    output_dir: str | Path,
    max_dim: int = 4000,
) -> dict:
    """Write one mask-based cavity QC image set without per-bin rendering.

    The binSegment 3-5 overlays use the continuous, nearest-neighbor-expanded
    GEM canvas. Reusing those images keeps zero-expression cells inside the GEM
    mask colored while the cavity tissue mask controls the final visible area.
    """
    source_overlay_dir = Path(source_overlay_dir)
    tissue_mask_path = Path(tissue_mask_path)
    output_dir = Path(output_dir)
    if not source_overlay_dir.is_dir():
        raise FileNotFoundError(f"Source overlay directory not found: {source_overlay_dir}")
    if not tissue_mask_path.is_file():
        raise FileNotFoundError(f"Cavity tissue mask not found: {tissue_mask_path}")

    source_files = sorted(
        path
        for path in source_overlay_dir.iterdir()
        if path.is_file()
        and path.name.startswith(_CAVITY_OVERLAY_PREFIXES)
        and ".full_canvas_qc." not in path.name
    )
    if not source_files:
        raise FileNotFoundError(f"No binSegment 3-5 overlay images found in: {source_overlay_dir}")

    with Image.open(tissue_mask_path) as mask_image:
        native_mask = mask_image.convert("L").copy()
    canvas_source = next(
        (path for path in source_files if path.name.startswith("3c_")),
        source_files[0],
    )
    with Image.open(canvas_source) as canvas_image:
        canvas_height, canvas_width = _resize_shape(
            (canvas_image.height, canvas_image.width),
            max_dim=max_dim,
        )
    canvas_size = (canvas_width, canvas_height)

    with Image.open(canvas_source) as heatmap_image:
        continuous_heatmap_pixels = np.asarray(heatmap_image.convert("RGB")).copy()
    continuous_heatmap_pixels[np.all(continuous_heatmap_pixels >= 250, axis=2)] = _JET_ZERO_RGB
    continuous_heatmap = Image.fromarray(continuous_heatmap_pixels).resize(
        canvas_size,
        resample=Image.Resampling.LANCZOS,
    )

    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    written_files: list[str] = []
    mask_output = output_dir / "3_tissue_segmentation_mask.png"
    display_mask = native_mask.resize(canvas_size, resample=Image.Resampling.NEAREST)
    display_mask.save(mask_output)
    written_files.append(mask_output.name)

    for source_path in source_files:
        destination = output_dir / source_path.name
        with Image.open(source_path) as source_image:
            if source_path.name.startswith("3a_"):
                # In the cavity output, 3a must show the effective cavity
                # mask instead of repeating binSegment's upstream HE mask.
                display_mask.save(destination)
                written_files.append(destination.name)
                continue
            if source_path.name.startswith("3c_"):
                rendered = continuous_heatmap.copy()
            elif source_path.name.startswith("5_"):
                background_name = "2_ssdna_registered.jpg" if "_on_ssdna" in source_path.name else "2_he_registered.jpg"
                background_path = source_overlay_dir / background_name
                if background_path.exists():
                    with Image.open(background_path) as background_image:
                        background = background_image.convert("RGB").resize(
                            canvas_size,
                            resample=Image.Resampling.LANCZOS,
                        )
                    rendered = Image.blend(background, continuous_heatmap, alpha=0.4)
                else:
                    rendered = source_image.convert("RGB").resize(
                        canvas_size,
                        resample=Image.Resampling.LANCZOS,
                    )
            else:
                rendered = source_image.convert("RGB").resize(
                    canvas_size,
                    resample=Image.Resampling.LANCZOS,
                )
        resized_mask = native_mask.resize(rendered.size, resample=Image.Resampling.NEAREST)
        white = Image.new("RGB", rendered.size, color=(255, 255, 255))
        rendered = Image.composite(rendered, white, resized_mask)
        save_options = {"quality": 95} if destination.suffix.lower() in {".jpg", ".jpeg"} else {}
        rendered.save(destination, **save_options)
        written_files.append(destination.name)

    return {
        "status": "completed",
        "output_dir": str(output_dir),
        "source_overlay_dir": str(source_overlay_dir),
        "tissue_mask": str(tissue_mask_path),
        "rendering": "binSegment continuous GEM overlays clipped by cavity tissue mask",
        "files": written_files,
        "max_dim": int(max_dim),
    }


def _coordinate_shape(positions: dict[str, tuple[float, float]]) -> tuple[int, int]:
    max_row = max((coord[0] for coord in positions.values()), default=1023.0)
    max_col = max((coord[1] for coord in positions.values()), default=1023.0)
    return max(1, int(np.ceil(max_row)) + 1), max(1, int(np.ceil(max_col)) + 1)


def _source_registered_shape(source_overlay_dir: Path | None) -> tuple[tuple[int, int] | None, str | None]:
    if source_overlay_dir is None:
        return None, None
    image_dir = source_overlay_dir.parent
    candidates = []
    for pattern in ("*_regist.tif", "*_regist.tiff", "*_regist.png"):
        candidates.extend(sorted(image_dir.glob(pattern)))
    for path in candidates:
        shape = _image_shape(path)
        if shape is not None:
            return shape, str(path)
    return None, None


def _shape_fit_score(coord_shape: tuple[int, int], image_shape: tuple[int, int]) -> float:
    return max(coord_shape[0] / max(1.0, float(image_shape[0])), coord_shape[1] / max(1.0, float(image_shape[1])))


def _resolve_shapes(
    output_bin_dir: Path,
    positions: dict[str, tuple[float, float]],
    source_overlay_dir: Path | None = None,
    max_dim: int = 4000,
) -> tuple[tuple[int, int], tuple[int, int], str]:
    hires_shape = _image_shape(output_bin_dir / "spatial" / "tissue_hires_image.png")
    coord_shape = _coordinate_shape(positions)
    registered_shape, registered_source = _source_registered_shape(source_overlay_dir)

    if hires_shape is not None and _shape_fit_score(coord_shape, hires_shape) <= 1.05:
        reference_shape = hires_shape
        reference_source = str(output_bin_dir / "spatial" / "tissue_hires_image.png")
    elif registered_shape is not None and _shape_fit_score(coord_shape, registered_shape) <= 1.05:
        reference_shape = registered_shape
        reference_source = registered_source or "registered_image"
    else:
        reference_shape = coord_shape
        reference_source = "position_extent"

    overlay_shape = None
    if source_overlay_dir is not None:
        overlay_shape = _image_shape(source_overlay_dir / "1_gem_expression.png")
    if overlay_shape is not None:
        target_shape = overlay_shape
        shape_source = f"{source_overlay_dir / '1_gem_expression.png'}; reference={reference_source}"
    else:
        target_shape = _resize_shape(reference_shape, max_dim=max_dim)
        shape_source = reference_source

    target_shape = _resize_shape(target_shape, max_dim=max_dim)
    return reference_shape, target_shape, shape_source


def _robust_normalize(values: np.ndarray) -> np.ndarray:
    out = np.zeros(values.shape, dtype=np.uint8)
    positive = values > 0
    if not positive.any():
        return out
    log_values = np.zeros_like(values, dtype=np.float32)
    log_values[positive] = np.log1p(values[positive])
    positive_values = log_values[positive]
    high = float(np.percentile(positive_values, 99.5))
    low = float(np.percentile(positive_values, 1.0))
    if high <= low:
        high = float(positive_values.max())
        low = float(positive_values.min())
    if high <= low:
        out[positive] = 255
        return out
    scaled = (log_values - low) / (high - low)
    return np.clip(scaled * 255.0, 0, 255).astype(np.uint8)


def _interpolate_colormap(gray_u8: np.ndarray, colors: np.ndarray) -> np.ndarray:
    gray = np.asarray(gray_u8, dtype=np.float32) / 255.0
    x = np.linspace(0.0, 1.0, colors.shape[0], dtype=np.float32)
    out = np.empty((*gray.shape, 3), dtype=np.uint8)
    for channel in range(3):
        out[:, :, channel] = np.interp(gray, x, colors[:, channel]).astype(np.uint8)
    return out


def _viridis(gray_u8: np.ndarray) -> np.ndarray:
    colors = np.array(
        [
            [68, 1, 84],
            [59, 82, 139],
            [33, 145, 140],
            [94, 201, 98],
            [253, 231, 37],
        ],
        dtype=np.float32,
    )
    return _interpolate_colormap(gray_u8, colors)


def _jet(gray_u8: np.ndarray) -> np.ndarray:
    colors = np.array(
        [
            [0, 0, 128],
            [0, 0, 255],
            [0, 255, 255],
            [255, 255, 0],
            [255, 0, 0],
            [128, 0, 0],
        ],
        dtype=np.float32,
    )
    return _interpolate_colormap(gray_u8, colors)


def _infer_bin_step(coords: np.ndarray) -> float:
    """Infer the native center-to-center spacing along one spatial axis."""
    unique = np.unique(coords)
    if unique.size < 2:
        return 1.0
    delta = np.diff(unique)
    delta = delta[delta > 0]
    return float(np.percentile(delta, 10)) if delta.size else 1.0


def _render_canvas(
    barcodes: list[str],
    umi_values: np.ndarray,
    positions: dict[str, tuple[float, float]],
    reference_shape: tuple[int, int],
    target_shape: tuple[int, int],
) -> np.ndarray:
    target_h, target_w = target_shape
    ref_h, ref_w = reference_shape
    canvas = np.zeros((target_h, target_w), dtype=np.float32)
    if not barcodes:
        return canvas

    coords = [positions[barcode] for barcode in barcodes if barcode in positions]
    values = [float(umi_values[idx]) for idx, barcode in enumerate(barcodes) if barcode in positions]
    if not coords:
        return canvas

    coord_array = np.asarray(coords, dtype=np.float32)
    values_array = np.asarray(values, dtype=np.float32)
    scale_y = target_h / max(1.0, float(ref_h))
    scale_x = target_w / max(1.0, float(ref_w))
    y = np.rint(coord_array[:, 0] * scale_y).astype(np.int32)
    x = np.rint(coord_array[:, 1] * scale_x).astype(np.int32)
    valid = (y >= 0) & (y < target_h) & (x >= 0) & (x < target_w)
    coord_array = coord_array[valid]
    values_array = values_array[valid]
    step_y = _infer_bin_step(coord_array[:, 0])
    step_x = _infer_bin_step(coord_array[:, 1])

    for (yy, xx), value in zip(coord_array, values_array):
        # Project each native bin's physical square instead of expanding a
        # center point. Adjacent bins meet at their shared edge and never
        # overlap; a genuinely absent bin remains visibly blank.
        y0 = max(0, int(np.rint((yy - step_y / 2.0) * scale_y)))
        y1 = min(target_h, int(np.rint((yy + step_y / 2.0) * scale_y)))
        x0 = max(0, int(np.rint((xx - step_x / 2.0) * scale_x)))
        x1 = min(target_w, int(np.rint((xx + step_x / 2.0) * scale_x)))
        if y1 <= y0:
            y1 = min(target_h, y0 + 1)
        if x1 <= x0:
            x1 = min(target_w, x0 + 1)
        canvas[y0:y1, x0:x1] = value
    return canvas


def write_filtered_gem_qc(
    output_bin_dir: str | Path,
    output_dir: str | Path,
    source_overlay_dir: str | Path | None = None,
    max_dim: int = 4000,
) -> dict:
    """Write filtered GEM-expression QC images for one square-bin output."""
    output_bin_dir = Path(output_bin_dir)
    output_dir = Path(output_dir)
    source_overlay_dir = Path(source_overlay_dir) if source_overlay_dir else None

    matrix_dir = output_bin_dir / "filtered_feature_bc_matrix"
    spatial_dir = output_bin_dir / "spatial"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    positions_path = spatial_dir / "tissue_positions_list.csv"
    for required in [barcodes_path, matrix_path, positions_path]:
        if not required.exists():
            raise FileNotFoundError(f"Required input not found: {required}")

    barcodes = _load_barcodes(barcodes_path)
    positions = _load_positions(positions_path)
    matrix = sio.mmread(str(matrix_path)).tocsr()
    umi_values = np.asarray(matrix.sum(axis=0)).ravel()

    reference_shape, target_shape, shape_source = _resolve_shapes(
        output_bin_dir=output_bin_dir,
        positions=positions,
        source_overlay_dir=source_overlay_dir,
        max_dim=int(max_dim),
    )
    canvas = _render_canvas(
        barcodes=barcodes,
        umi_values=umi_values,
        positions=positions,
        reference_shape=reference_shape,
        target_shape=target_shape,
    )
    normalized = _robust_normalize(canvas)
    expression_rgb = _viridis(normalized)

    positive = canvas > 0
    heatmap_rgb = np.full((*normalized.shape, 3), 255, dtype=np.uint8)
    jet_rgb = _jet(normalized)
    heatmap_rgb[positive] = jet_rgb[positive]

    output_dir.mkdir(parents=True, exist_ok=True)
    expression_path = output_dir / "1_filtered_gem_expression.png"
    heatmap_path = output_dir / "2_filtered_gem_heatmap_only.png"
    Image.fromarray(expression_rgb).save(expression_path)
    Image.fromarray(heatmap_rgb).save(heatmap_path)

    summary = {
        "status": "completed",
        "output_dir": str(output_dir),
        "gem_expression": str(expression_path),
        "gem_heatmap_only": str(heatmap_path),
        "barcodes": int(len(barcodes)),
        "nonzero_pixels": int(positive.sum()),
        "reference_shape": list(reference_shape),
        "target_shape": list(target_shape),
        "shape_source": shape_source,
        "max_dim": int(max_dim),
    }
    with open(output_dir / "qc_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    return summary
