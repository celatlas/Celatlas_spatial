#!/usr/bin/env python3
"""Run StarDist nucleus segmentation and write labels/QC overlays."""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

Image.MAX_IMAGE_PIXELS = None


FLUORESCENCE_NAME_HINTS = ("dapi", "fluores", "fluoro", "immunofluores", "multiplex", "mif")


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run StarDist 2D nucleus segmentation.")
    parser.add_argument("--image", required=True, help="Registered tissue image (HE or fluorescence).")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--sample", required=True, help="Sample name used in output filenames.")
    parser.add_argument(
        "--model",
        default=None,
        help=(
            "StarDist pretrained model name. Defaults to 2D_versatile_fluo for fluorescence-like "
            "inputs and 2D_versatile_he otherwise."
        ),
    )
    parser.add_argument(
        "--model-dir",
        help=(
            "Local StarDist model directory containing config.json, thresholds.json, and weights. "
            "If set, bypasses StarDist2D.from_pretrained."
        ),
    )
    parser.add_argument("--prob-thresh", type=float, default=0.30, help="StarDist probability threshold.")
    parser.add_argument("--nms-thresh", type=float, default=None, help="Optional StarDist NMS threshold.")
    parser.add_argument("--max-dim", type=int, default=6000, help="Resize image long edge to this size; 0 disables resizing.")
    parser.add_argument(
        "--scale",
        type=float,
        default=2.0,
        help="StarDist internal spatial scale. Use 2.0 to match the validated H&E workflow; use 1.0 to disable.",
    )
    parser.add_argument("--n-tiles", default="4,4", help="Comma-separated StarDist tile grid, e.g. 4,4.")
    parser.add_argument("--preview-max-dim", type=int, default=3000, help="Resize overlay preview long edge to this size.")
    parser.add_argument("--normalization-low", type=float, default=1.0)
    parser.add_argument("--normalization-high", type=float, default=99.8)
    return parser.parse_args()


def _default_model_name(image_path: Path | None) -> str:
    if image_path is not None:
        path_text = " ".join(part.lower() for part in image_path.parts)
        if any(hint in path_text for hint in FLUORESCENCE_NAME_HINTS):
            return "2D_versatile_fluo"
    return "2D_versatile_he"


def _default_model_dir(model_name: str | None = None) -> str | None:
    for env_name in ("STARDIST_MODEL_DIR", "CELATLAS_STARDIST_MODEL_DIR"):
        env_value = os.environ.get(env_name)
        if env_value:
            return env_value

    candidate_models = [model_name] if model_name else ["2D_versatile_he", "2D_versatile_fluo"]
    for candidate in candidate_models:
        if not candidate:
            continue
        local_cache = (
            Path.home()
            / ".keras"
            / "models"
            / "StarDist2D"
            / candidate
            / f"{candidate}_extracted"
        )
        if local_cache.exists():
            return str(local_cache)
    return None


def resize_image(image: Image.Image, max_dim: int) -> tuple[Image.Image, float]:
    if max_dim <= 0 or max(image.size) <= max_dim:
        return image, 1.0
    scale = float(max_dim) / float(max(image.size))
    new_size = (max(1, int(round(image.size[0] * scale))), max(1, int(round(image.size[1] * scale))))
    try:
        resample = Image.Resampling.BILINEAR
    except AttributeError:
        resample = Image.BILINEAR
    return image.resize(new_size, resample), scale


def parse_tiles(value: str, ndim: int) -> tuple[int, ...]:
    parts = [int(part.strip()) for part in value.split(",") if part.strip()]
    if len(parts) == 2 and ndim == 3:
        parts.append(1)
    if len(parts) != ndim or any(part < 1 for part in parts):
        raise ValueError(f"Invalid --n-tiles value: {value}")
    return tuple(parts)


def _as_preview_rgb(image: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        return np.repeat(image[..., None], 3, axis=2)
    if image.ndim == 3 and image.shape[2] == 1:
        return np.repeat(image, 3, axis=2)
    if image.ndim == 3 and image.shape[2] >= 3:
        return image[..., :3]
    raise ValueError(f"Unsupported image shape for preview: {image.shape}")


def _prepare_image_for_model(image: Image.Image, n_channel_in: int) -> Image.Image:
    if int(n_channel_in) <= 1:
        return image.convert("L")
    return image.convert("RGB")


def normalize_model_input(
    image: np.ndarray,
    pmin: float,
    pmax: float,
    *,
    background_delta: float = 2.0,
) -> tuple[np.ndarray, dict]:
    """Normalize an image without letting a sparse uniform canvas collapse its range."""
    image = np.asarray(image)
    spatial_axes = tuple(range(min(2, image.ndim)))
    low = np.percentile(image, pmin, axis=spatial_axes, keepdims=True)
    high = np.percentile(image, pmax, axis=spatial_axes, keepdims=True)
    span = high - low

    summary = {
        "method": "full_image_percentile",
        "low": np.asarray(low).reshape(-1).astype(float).tolist(),
        "high": np.asarray(high).reshape(-1).astype(float).tolist(),
        "foreground_fraction": None,
    }
    if np.all(span >= 1.0):
        normalized = (image.astype(np.float32) - low.astype(np.float32)) / span.astype(np.float32)
        return normalized, summary

    # A registered tissue image can contain less than pmin percent tissue. In
    # that case both percentiles equal the white (or black) canvas value and
    # normalization amplifies tissue pixels into enormous model activations.
    if image.ndim == 2:
        background = float(np.median(image))
        foreground = np.abs(image.astype(np.float32) - background) > background_delta
        foreground_values = image[foreground][:, None]
    else:
        background = np.median(image, axis=spatial_axes, keepdims=True)
        foreground = np.zeros(image.shape[:2], dtype=bool)
        for channel in range(image.shape[2]):
            difference = np.abs(image[..., channel].astype(np.float32) - float(background[..., channel].item()))
            foreground |= difference > background_delta
        foreground_values = image[foreground]

    foreground_fraction = float(np.mean(foreground))
    if foreground_values.shape[0] < 64:
        raise ValueError(
            "StarDist input has no usable intensity range: "
            f"full-image percentile span={np.asarray(span).reshape(-1).tolist()}, "
            f"foreground_pixels={foreground_values.shape[0]}"
        )

    fallback_low = np.percentile(foreground_values, pmin, axis=0)
    fallback_high = np.percentile(foreground_values, pmax, axis=0)
    fallback_span = fallback_high - fallback_low
    if np.any(fallback_span < 1.0):
        fallback_low = np.min(foreground_values, axis=0)
        fallback_high = np.max(foreground_values, axis=0)
        fallback_span = fallback_high - fallback_low
    if np.any(fallback_span < 1.0):
        raise ValueError(
            "StarDist foreground has no usable intensity range: "
            f"range={np.asarray(fallback_span).reshape(-1).tolist()}"
        )

    if image.ndim == 2:
        fallback_low = float(np.asarray(fallback_low).reshape(-1)[0])
        fallback_span = float(np.asarray(fallback_span).reshape(-1)[0])
    else:
        fallback_low = np.asarray(fallback_low).reshape((1, 1, -1))
        fallback_span = np.asarray(fallback_span).reshape((1, 1, -1))
    normalized = (image.astype(np.float32) - np.asarray(fallback_low, dtype=np.float32)) / np.asarray(
        fallback_span, dtype=np.float32
    )
    normalized = np.clip(normalized, 0.0, 1.0)
    summary.update(
        {
            "method": "foreground_percentile",
            "low": np.asarray(fallback_low).reshape(-1).astype(float).tolist(),
            "high": np.asarray(fallback_high).reshape(-1).astype(float).tolist(),
            "foreground_fraction": foreground_fraction,
        }
    )
    return normalized, summary


def write_boundary_preview(image_rgb: np.ndarray, labels: np.ndarray, output_path: Path, max_dim: int) -> None:
    from skimage.segmentation import find_boundaries

    preview = _as_preview_rgb(image_rgb).copy()
    boundaries = find_boundaries(labels, mode="outer")
    preview[boundaries] = np.array([255, 40, 40], dtype=np.uint8)
    preview_image = Image.fromarray(preview)
    if max_dim > 0 and max(preview_image.size) > max_dim:
        scale = float(max_dim) / float(max(preview_image.size))
        new_size = (
            max(1, int(round(preview_image.size[0] * scale))),
            max(1, int(round(preview_image.size[1] * scale))),
        )
        try:
            resample = Image.Resampling.BILINEAR
        except AttributeError:
            resample = Image.BILINEAR
        preview_image = preview_image.resize(new_size, resample)
    preview_image.save(output_path, quality=92)


def main() -> None:
    args = parse_args()
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

    from stardist.models import StarDist2D
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    image_path = Path(args.image)
    selected_model = args.model or _default_model_name(image_path)
    model_dir_path: Path | None = None

    model_dir_value = args.model_dir or _default_model_dir(selected_model)
    if model_dir_value:
        model_dir_path = Path(model_dir_value).expanduser().resolve()
        log(f"Loading StarDist model from local directory: {model_dir_path}")
        for required_name in ("config.json", "thresholds.json", "weights_best.h5"):
            required_path = model_dir_path / required_name
            if not required_path.exists():
                raise FileNotFoundError(f"Local StarDist model is missing {required_name}: {required_path}")
        model = StarDist2D(None, name=model_dir_path.name, basedir=str(model_dir_path.parent))
        model_name = model_dir_path.name
    else:
        log(f"Loading StarDist model: {selected_model}")
        model = StarDist2D.from_pretrained(selected_model)
        model_name = selected_model

    input_channels = int(getattr(model.config, "n_channel_in", 1) or 1)

    log(f"Loading image: {image_path}")
    image = Image.open(image_path)
    original_size = image.size
    image = _prepare_image_for_model(image, input_channels)
    image, resize_scale = resize_image(image, int(args.max_dim))
    image_array = np.asarray(image, dtype=np.uint8)
    preview_rgb = _as_preview_rgb(image_array)
    log(
        f"image_size={image.size}, original_size={original_size}, "
        f"resize_scale={resize_scale:.6f}, stardist_scale={args.scale:g}, "
        f"input_channels={input_channels}"
    )

    normalized, normalization = normalize_model_input(
        image_array,
        pmin=float(args.normalization_low),
        pmax=float(args.normalization_high),
    )
    log(
        f"normalization={normalization['method']}, low={normalization['low']}, "
        f"high={normalization['high']}, foreground_fraction={normalization['foreground_fraction']}"
    )
    n_tiles = parse_tiles(args.n_tiles, normalized.ndim)

    predict_kwargs = {"prob_thresh": float(args.prob_thresh), "n_tiles": n_tiles}
    if args.nms_thresh is not None:
        predict_kwargs["nms_thresh"] = float(args.nms_thresh)
    if float(args.scale) != 1.0:
        predict_kwargs["scale"] = float(args.scale)

    start = time.time()
    log(f"Running StarDist: prob_thresh={args.prob_thresh}, scale={args.scale:g}, n_tiles={n_tiles}")
    labels, _ = model.predict_instances(normalized, **predict_kwargs)
    elapsed_min = (time.time() - start) / 60.0
    labels = labels.astype(np.uint32, copy=False)
    label_count = int(labels.max())
    log(f"labels={label_count:,}, elapsed={elapsed_min:.1f} min")

    safe_model_name = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in model_name)
    scale_tag = f"scale{args.scale:g}".replace(".", "p")
    prefix = (
        f"{args.sample}.stardist_{safe_model_name}.maxdim{args.max_dim}_"
        f"{scale_tag}_p{int(round(args.prob_thresh * 100)):02d}"
    )
    labels_path = outdir / f"{prefix}.labels.tif"
    overlay_path = outdir / f"{prefix}.boundary_overlay.preview.jpg"
    summary_path = outdir / f"{prefix}.summary.json"
    bbox_path = outdir / f"{args.sample}.stardist_bbox.csv"
    latest_labels = outdir / f"{args.sample}.stardist_labels.tif"
    latest_overlay = outdir / f"{args.sample}.stardist_boundary_overlay.preview.jpg"

    tifffile.imwrite(labels_path, labels)
    tifffile.imwrite(latest_labels, labels)
    bbox_path.write_text(f"0,0,{original_size[0]},{original_size[1]}\n")
    write_boundary_preview(preview_rgb, labels, overlay_path, int(args.preview_max_dim))
    write_boundary_preview(preview_rgb, labels, latest_overlay, int(args.preview_max_dim))

    summary = {
        "sample": args.sample,
        "image": str(image_path.resolve()),
        "model": model_name,
        "model_dir": str(model_dir_path) if model_dir_path else None,
        "input_channels": input_channels,
        "input_image_mode": image.mode,
        "normalization": normalization,
        "prob_thresh": float(args.prob_thresh),
        "nms_thresh": args.nms_thresh,
        "n_tiles": list(n_tiles),
        "max_dim": int(args.max_dim),
        "resize_scale": float(resize_scale),
        "scale": float(args.scale),
        "stardist_scale": float(args.scale),
        "original_size": list(original_size),
        "processed_size": list(image.size),
        "labels": label_count,
        "labels_path": str(labels_path),
        "latest_labels": str(latest_labels),
        "bbox": str(bbox_path),
        "boundary_overlay_preview": str(overlay_path),
        "elapsed_min": elapsed_min,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
    (outdir / f"{args.sample}.stardist_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    log(f"labels_path={labels_path}")
    log(f"summary={summary_path}")


if __name__ == "__main__":
    main()
