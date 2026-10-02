"""Generate smooth binary tissue/cavity masks from registered tissue images."""

from __future__ import annotations

from pathlib import Path

import cv2
import numpy as np
from PIL import Image
from scipy import ndimage as ndi

Image.MAX_IMAGE_PIXELS = None

try:
    from skimage import morphology as sk_morphology
except Exception:  # pragma: no cover - optional dependency fallback
    sk_morphology = None


MASK_PRESETS = ("default", "adipose")
ADIPOSE_HOLE_AREA_FRACTION = 0.006
ADIPOSE_BLANK_CAVITY_AREA_FRACTION = 0.0002
ADIPOSE_BLANK_CAVITY_MIN_AREA = 50000
ADIPOSE_BLANK_CAVITY_GRAY_MIN = 245
ADIPOSE_BLANK_CAVITY_SATURATION_MAX = 20
ADIPOSE_BLANK_CAVITY_FRACTION = 0.98
ADIPOSE_BLANK_CAVITY_INNER_MARGIN = 25


def _ensure_odd(value: int) -> int:
    return value if value % 2 == 1 else value + 1


def _otsu_threshold_u8(image_u8: np.ndarray) -> float:
    threshold, _ = cv2.threshold(np.asarray(image_u8, dtype=np.uint8), 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return float(threshold)


def _remove_small_objects(mask: np.ndarray, min_size: int) -> np.ndarray:
    labeled, num = ndi.label(mask)
    if num == 0:
        return mask
    sizes = ndi.sum(mask, labeled, index=np.arange(1, num + 1))
    keep = np.zeros(num + 1, dtype=bool)
    keep[1:] = sizes >= min_size
    return keep[labeled]


def _remove_small_holes(mask: np.ndarray, area_threshold: int, method: str = "scipy") -> np.ndarray:
    if area_threshold <= 0:
        return mask
    if method == "skimage" and sk_morphology is not None:
        return sk_morphology.remove_small_holes(mask.astype(bool), area_threshold=int(area_threshold), connectivity=1)

    inv = ~mask
    labeled, num = ndi.label(inv)
    if num == 0:
        return mask
    sizes = ndi.sum(inv, labeled, index=np.arange(1, num + 1))
    out = mask.copy()
    border_labels = set(np.unique(np.concatenate([labeled[0, :], labeled[-1, :], labeled[:, 0], labeled[:, -1]])))
    for label_id, size in enumerate(sizes, start=1):
        if label_id in border_labels:
            continue
        if size < area_threshold:
            out[labeled == label_id] = True
    return out


def _keep_significant_components(mask: np.ndarray, min_ratio: float) -> np.ndarray:
    labeled, num = ndi.label(mask)
    if num == 0:
        return mask
    sizes = ndi.sum(mask, labeled, index=np.arange(1, num + 1))
    max_size = float(np.max(sizes))
    keep = np.zeros(num + 1, dtype=bool)
    keep[1:] = sizes >= max(0.0, float(min_ratio)) * max_size
    return keep[labeled]


def _find_blank_internal_holes(
    mask: np.ndarray,
    gray_u8: np.ndarray,
    saturation_u8: np.ndarray,
    min_area: int,
    gray_min: int,
    saturation_max: int,
    blank_fraction: float,
) -> np.ndarray:
    """Return holes that are enclosed by tissue and are nearly pure blank slide."""
    if min_area <= 0:
        return np.zeros(mask.shape, dtype=bool)

    inv = ~mask
    labeled, num = ndi.label(inv)
    if num == 0:
        return np.zeros(mask.shape, dtype=bool)

    index = np.arange(1, num + 1)
    sizes = np.asarray(ndi.sum(inv, labeled, index=index), dtype=float)
    blank_pixels = ((gray_u8 >= int(gray_min)) & (saturation_u8 <= int(saturation_max)) & inv).astype(np.uint8)
    blank_counts = np.asarray(ndi.sum(blank_pixels, labeled, index=index), dtype=float)
    fractions = np.divide(blank_counts, sizes, out=np.zeros_like(blank_counts), where=sizes > 0)

    border_labels = set(np.unique(np.concatenate([labeled[0, :], labeled[-1, :], labeled[:, 0], labeled[:, -1]])))
    keep_labels = [
        label_id
        for label_id, size, fraction in zip(index, sizes, fractions)
        if label_id not in border_labels and size >= int(min_area) and fraction >= float(blank_fraction)
    ]
    if not keep_labels:
        return np.zeros(mask.shape, dtype=bool)
    return np.isin(labeled, np.asarray(keep_labels, dtype=labeled.dtype))


def _find_large_blank_regions(
    mask: np.ndarray,
    gray_u8: np.ndarray,
    saturation_u8: np.ndarray,
    min_area: int,
    gray_min: int,
    saturation_max: int,
    inner_margin: int = 0,
) -> np.ndarray:
    """Return large blank-looking connected regions inside the current tissue mask."""
    if min_area <= 0:
        return np.zeros(mask.shape, dtype=bool)

    inner_mask = mask
    if inner_margin > 0:
        radius = int(inner_margin)
        kernel = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE,
            (_ensure_odd(radius * 2 + 1), _ensure_odd(radius * 2 + 1)),
        )
        inner_mask = cv2.erode(mask.astype(np.uint8) * 255, kernel) > 0

    blank = (gray_u8 >= int(gray_min)) & (saturation_u8 <= int(saturation_max)) & inner_mask
    labeled, num = ndi.label(blank)
    if num == 0:
        return np.zeros(mask.shape, dtype=bool)

    index = np.arange(1, num + 1)
    sizes = np.asarray(ndi.sum(blank, labeled, index=index), dtype=float)
    border_labels = set(np.unique(np.concatenate([labeled[0, :], labeled[-1, :], labeled[:, 0], labeled[:, -1]])))
    keep_labels = [
        label_id
        for label_id, size in zip(index, sizes)
        if label_id not in border_labels and size >= int(min_area)
    ]
    if not keep_labels:
        return np.zeros(mask.shape, dtype=bool)
    return np.isin(labeled, np.asarray(keep_labels, dtype=labeled.dtype))


def _resolve_mask_parameters(
    image_shape: tuple[int, ...],
    blur_kernel: int,
    close_radius: int,
    open_radius: int,
    min_object_area: int,
    min_hole_area: int,
    min_component_ratio: float,
    mask_preset: str = "default",
) -> tuple[dict, dict]:
    mask_preset = (mask_preset or "default").lower()
    if mask_preset not in MASK_PRESETS:
        raise ValueError(f"Unsupported mask preset: {mask_preset}")

    params = {
        "blur_kernel": int(blur_kernel),
        "close_radius": int(close_radius),
        "open_radius": int(open_radius),
        "min_object_area": int(min_object_area),
        "min_hole_area": int(min_hole_area),
        "min_component_ratio": float(min_component_ratio),
        "hole_fill_method": "scipy",
        "blank_cavity_min_area": 0,
        "blank_cavity_gray_min": ADIPOSE_BLANK_CAVITY_GRAY_MIN,
        "blank_cavity_saturation_max": ADIPOSE_BLANK_CAVITY_SATURATION_MAX,
        "blank_cavity_fraction": ADIPOSE_BLANK_CAVITY_FRACTION,
        "blank_cavity_inner_margin": 0,
    }
    notes = {
        "mask_preset": mask_preset,
        "requested_blur_kernel": int(blur_kernel),
        "requested_close_radius": int(close_radius),
        "requested_open_radius": int(open_radius),
        "requested_min_object_area": int(min_object_area),
        "requested_min_hole_area": int(min_hole_area),
        "requested_min_component_ratio": float(min_component_ratio),
    }

    if mask_preset == "adipose":
        height, width = int(image_shape[0]), int(image_shape[1])
        adaptive_hole_area = max(1, int(round(height * width * ADIPOSE_HOLE_AREA_FRACTION)))
        params["open_radius"] = 0
        params["min_component_ratio"] = 0.0
        params["min_hole_area"] = max(params["min_hole_area"], adaptive_hole_area)
        params["hole_fill_method"] = "skimage" if sk_morphology is not None else "scipy"
        params["blank_cavity_min_area"] = max(
            ADIPOSE_BLANK_CAVITY_MIN_AREA,
            int(round(height * width * ADIPOSE_BLANK_CAVITY_AREA_FRACTION)),
        )
        params["blank_cavity_inner_margin"] = ADIPOSE_BLANK_CAVITY_INNER_MARGIN
        notes["adipose_hole_area_fraction"] = float(ADIPOSE_HOLE_AREA_FRACTION)
        notes["adipose_adaptive_min_hole_area"] = int(adaptive_hole_area)
        notes["adipose_blank_cavity_area_fraction"] = float(ADIPOSE_BLANK_CAVITY_AREA_FRACTION)
        notes["skimage_available"] = sk_morphology is not None

    return params, notes


def build_tissue_mask(
    image_rgb: np.ndarray,
    blur_kernel: int = 5,
    close_radius: int = 7,
    open_radius: int = 1,
    min_object_area: int = 3000,
    min_hole_area: int = 1500,
    min_component_ratio: float = 0.0,
    mask_preset: str = "default",
) -> np.ndarray:
    """Return a boolean mask where True means valid tissue."""
    image_rgb = np.asarray(image_rgb, dtype=np.uint8)
    params, _ = _resolve_mask_parameters(
        image_shape=image_rgb.shape,
        blur_kernel=blur_kernel,
        close_radius=close_radius,
        open_radius=open_radius,
        min_object_area=min_object_area,
        min_hole_area=min_hole_area,
        min_component_ratio=min_component_ratio,
        mask_preset=mask_preset,
    )
    blur_kernel = params["blur_kernel"]
    close_radius = params["close_radius"]
    open_radius = params["open_radius"]
    min_object_area = params["min_object_area"]
    min_hole_area = params["min_hole_area"]
    min_component_ratio = params["min_component_ratio"]
    hole_fill_method = params["hole_fill_method"]
    blank_cavity_min_area = params["blank_cavity_min_area"]

    hsv = cv2.cvtColor(image_rgb, cv2.COLOR_RGB2HSV)
    gray = cv2.cvtColor(image_rgb, cv2.COLOR_RGB2GRAY)

    blur_kernel = _ensure_odd(max(3, int(blur_kernel)))
    gray_blur = cv2.GaussianBlur(gray, (blur_kernel, blur_kernel), 0)

    saturation_u8 = hsv[:, :, 1]
    sat_thresh = _otsu_threshold_u8(saturation_u8)
    val_thresh = _otsu_threshold_u8(gray_blur)
    tissue = (saturation_u8 > max(8, int(sat_thresh * 0.65))) | (gray_blur < min(245, int(val_thresh * 1.05)))

    kernel_open = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE,
        (_ensure_odd(max(1, int(open_radius) * 2 + 1)), _ensure_odd(max(1, int(open_radius) * 2 + 1))),
    )
    kernel_close = cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE,
        (_ensure_odd(max(1, int(close_radius) * 2 + 1)), _ensure_odd(max(1, int(close_radius) * 2 + 1))),
    )

    if open_radius > 0:
        tissue = cv2.morphologyEx(tissue.astype(np.uint8) * 255, cv2.MORPH_OPEN, kernel_open) > 0
    if close_radius > 0:
        tissue = cv2.morphologyEx(tissue.astype(np.uint8) * 255, cv2.MORPH_CLOSE, kernel_close) > 0

    tissue = _remove_small_objects(tissue, int(min_object_area))
    tissue = _keep_significant_components(tissue, float(min_component_ratio))
    blank_cavity = _find_blank_internal_holes(
        tissue,
        gray_u8=gray,
        saturation_u8=saturation_u8,
        min_area=int(blank_cavity_min_area),
        gray_min=int(params["blank_cavity_gray_min"]),
        saturation_max=int(params["blank_cavity_saturation_max"]),
        blank_fraction=float(params["blank_cavity_fraction"]),
    )
    tissue = _remove_small_holes(tissue, int(min_hole_area), method=hole_fill_method)
    if blank_cavity.any():
        tissue[blank_cavity] = False
    blank_regions = _find_large_blank_regions(
        tissue,
        gray_u8=gray,
        saturation_u8=saturation_u8,
        min_area=int(blank_cavity_min_area),
        gray_min=int(params["blank_cavity_gray_min"]),
        saturation_max=int(params["blank_cavity_saturation_max"]),
        inner_margin=int(params["blank_cavity_inner_margin"]),
    )
    if blank_regions.any():
        tissue[blank_regions] = False
    return tissue


def build_cavity_mask(tissue_mask: np.ndarray, min_object_area: int = 3000) -> np.ndarray:
    """Return a boolean mask where True means inner cavity/hole."""
    filled = ndi.binary_fill_holes(tissue_mask)
    cavity = filled & (~tissue_mask)
    return _remove_small_objects(cavity, max(64, int(min_object_area) // 4))


def generate_binary_mask(
    input_image: str | Path,
    output_mask: str | Path,
    mode: str = "tissue",
    max_dim: int | None = 4096,
    blur_kernel: int = 5,
    close_radius: int = 7,
    open_radius: int = 1,
    min_object_area: int = 3000,
    min_hole_area: int = 1500,
    min_component_ratio: float = 0.0,
    mask_preset: str = "default",
) -> dict:
    """Generate and write a binary mask. White pixels are foreground."""
    input_image = Path(input_image)
    output_mask = Path(output_mask)
    try:
        image = Image.open(input_image).convert("RGB")
    except Exception as exc:
        raise FileNotFoundError(f"Failed to read image: {input_image}") from exc
    original_shape = (image.size[1], image.size[0])
    if max_dim and max_dim > 0 and max(image.size) > max_dim:
        scale = float(max_dim) / float(max(image.size))
        new_size = (max(1, int(round(image.size[0] * scale))), max(1, int(round(image.size[1] * scale))))
        try:
            resample = Image.Resampling.BILINEAR
        except AttributeError:
            resample = Image.BILINEAR
        image = image.resize(new_size, resample)
    image_rgb = np.asarray(image, dtype=np.uint8)

    effective_params, preset_notes = _resolve_mask_parameters(
        image_shape=image_rgb.shape,
        blur_kernel=blur_kernel,
        close_radius=close_radius,
        open_radius=open_radius,
        min_object_area=min_object_area,
        min_hole_area=min_hole_area,
        min_component_ratio=min_component_ratio,
        mask_preset=mask_preset,
    )
    tissue_mask = build_tissue_mask(
        image_rgb=image_rgb,
        mask_preset=mask_preset,
        **{key: effective_params[key] for key in (
            "blur_kernel",
            "close_radius",
            "open_radius",
            "min_object_area",
            "min_hole_area",
            "min_component_ratio",
        )},
    )
    if mode == "tissue":
        mask = tissue_mask
    elif mode == "cavity":
        mask = build_cavity_mask(tissue_mask, min_object_area=min_object_area)
    else:
        raise ValueError(f"Unsupported mask mode: {mode}")

    output_mask.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(mask.astype(np.uint8) * 255).save(output_mask)
    return {
        "input_image": str(input_image),
        "output_mask": str(output_mask),
        "mode": mode,
        "shape": list(mask.shape),
        "original_image_shape": list(original_shape),
        "processed_image_shape": list(mask.shape),
        "max_dim": None if max_dim is None else int(max_dim),
        "foreground_pixels": int(mask.sum()),
        "foreground_fraction": float(mask.sum() / mask.size) if mask.size else 0.0,
        "mask_preset": preset_notes["mask_preset"],
        "blur_kernel": int(effective_params["blur_kernel"]),
        "close_radius": int(effective_params["close_radius"]),
        "open_radius": int(effective_params["open_radius"]),
        "min_object_area": int(effective_params["min_object_area"]),
        "min_hole_area": int(effective_params["min_hole_area"]),
        "min_component_ratio": float(effective_params["min_component_ratio"]),
        "hole_fill_method": effective_params["hole_fill_method"],
        "blank_cavity_min_area": int(effective_params["blank_cavity_min_area"]),
        "blank_cavity_gray_min": int(effective_params["blank_cavity_gray_min"]),
        "blank_cavity_saturation_max": int(effective_params["blank_cavity_saturation_max"]),
        "blank_cavity_fraction": float(effective_params["blank_cavity_fraction"]),
        "blank_cavity_inner_margin": int(effective_params["blank_cavity_inner_margin"]),
        **preset_notes,
    }
