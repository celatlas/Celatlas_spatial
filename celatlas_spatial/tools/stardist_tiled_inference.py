#!/usr/bin/env python3
"""Run StarDist on overlapping image tiles and rebuild a full-resolution label TIFF."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterator

import numpy as np
import tifffile
from PIL import Image

from celatlas_spatial.tools.stardist_inference import (
    _default_model_dir,
    _default_model_name,
    parse_tiles,
)


Image.MAX_IMAGE_PIXELS = None

NORMALIZATION_POLICY = "percentile-range-fallback-v1"


class ImageRegionSource:
    """Read rectangular regions without requiring callers to load the full image."""

    def __init__(self, image_path: Path):
        self.image_path = image_path
        self._array = self._open_contiguous_tiff_memmap(image_path)
        self._image: Image.Image | None = None
        if self._array is not None:
            self.height, self.width = self._array.shape[:2]
            self.backend = "tiff_memmap"
        else:
            self._image = Image.open(image_path)
            self.width, self.height = self._image.size
            self.backend = "pillow"

    @staticmethod
    def _open_contiguous_tiff_memmap(image_path: Path) -> np.memmap | None:
        if image_path.suffix.lower() not in {".tif", ".tiff"}:
            return None
        try:
            with tifffile.TiffFile(image_path) as tif:
                if len(tif.series) != 1 or len(tif.pages) != 1:
                    return None
                page = tif.pages[0]
                if page.is_tiled or int(page.compression) != 1 or int(page.predictor) != 1:
                    return None
                if page.axes not in {"YX", "YXS"}:
                    return None
                if page.planarconfig is not None and int(page.planarconfig) != 1:
                    return None
                orientation = page.tags.get("Orientation")
                if orientation is not None and int(orientation.value) != 1:
                    return None
                dtype = np.dtype(page.dtype)
                if dtype.itemsize > 1 and not dtype.isnative:
                    return None
                bits_per_sample = page.bitspersample
                if isinstance(bits_per_sample, tuple):
                    if any(int(bits) != dtype.itemsize * 8 for bits in bits_per_sample):
                        return None
                elif int(bits_per_sample) != dtype.itemsize * 8:
                    return None
                offsets = np.asarray(page.dataoffsets, dtype=np.int64)
                bytecounts = np.asarray(page.databytecounts, dtype=np.int64)
                if offsets.size == 0 or np.any(offsets[1:] != offsets[:-1] + bytecounts[:-1]):
                    return None
                height, width = page.shape[:2]
                samples = page.shape[2] if page.axes == "YXS" else 1
                row_bytes = width * samples * dtype.itemsize
                total_bytes = int(bytecounts.sum())
                if row_bytes <= 0 or total_bytes % row_bytes != 0:
                    return None
                stored_rows = total_bytes // row_bytes
                if stored_rows < height:
                    return None
                stored_shape = (stored_rows, width, samples) if samples > 1 else (stored_rows, width)
                data_offset = int(offsets[0])
            mapped = np.memmap(image_path, dtype=dtype, mode="r", offset=data_offset, shape=stored_shape)
            return mapped[:height]
        except (OSError, ValueError, tifffile.TiffFileError):
            return None

    @property
    def size(self) -> tuple[int, int]:
        return self.width, self.height

    def crop(self, box: tuple[int, int, int, int]) -> Image.Image:
        if self._array is not None:
            x0, y0, x1, y1 = box
            return Image.fromarray(np.array(self._array[y0:y1, x0:x1], copy=True))
        assert self._image is not None
        return self._image.crop(box)

    def preview_rgb(self, max_dim: int) -> Image.Image:
        return self.preview_crop_rgb((0, 0, self.width, self.height), max_dim)

    def preview_crop_rgb(
        self,
        box: tuple[int, int, int, int],
        max_dim: int,
    ) -> Image.Image:
        """Read an RGB preview of a region without materializing a whole slide."""
        x0, y0, x1, y1 = box
        if not (0 <= x0 < x1 <= self.width and 0 <= y0 < y1 <= self.height):
            raise ValueError(f"Invalid image crop: {box}")
        crop_width = x1 - x0
        crop_height = y1 - y0
        scale = min(1.0, float(max_dim) / max(crop_width, crop_height)) if max_dim > 0 else 1.0
        preview_size = (
            max(1, round(crop_width * scale)),
            max(1, round(crop_height * scale)),
        )
        if self._array is not None:
            y_index = y0 + np.minimum(
                (np.arange(preview_size[1], dtype=np.float64) / scale).astype(np.int64),
                crop_height - 1,
            )
            x_index = x0 + np.minimum(
                (np.arange(preview_size[0], dtype=np.float64) / scale).astype(np.int64),
                crop_width - 1,
            )
            sampled = np.asarray(self._array[np.ix_(y_index, x_index)])
            return Image.fromarray(sampled).convert("RGB")
        assert self._image is not None
        image = self._image.crop(box).convert("RGB")
        if image.size != preview_size:
            image = image.resize(preview_size, Image.Resampling.BILINEAR)
        return image

    def close(self) -> None:
        if self._image is not None:
            self._image.close()
        self._array = None

    def __enter__(self) -> "ImageRegionSource":
        return self

    def __exit__(self, *_args) -> None:
        self.close()


@dataclass(frozen=True)
class Tile:
    index: int
    row: int
    column: int
    core_x0: int
    core_y0: int
    core_x1: int
    core_y1: int
    read_x0: int
    read_y0: int
    read_x1: int
    read_y1: int

    @property
    def read_box(self) -> tuple[int, int, int, int]:
        return self.read_x0, self.read_y0, self.read_x1, self.read_y1

    @property
    def read_shape(self) -> tuple[int, int]:
        return self.read_y1 - self.read_y0, self.read_x1 - self.read_x0


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def iter_tiles(width: int, height: int, tile_size: int, overlap: int) -> Iterator[Tile]:
    """Yield non-overlapping cores with an overlapping inference halo."""
    if width <= 0 or height <= 0:
        raise ValueError("Image width and height must be positive")
    if tile_size <= 0:
        raise ValueError("--tile-size must be positive")
    if overlap < 0 or overlap * 2 >= tile_size:
        raise ValueError("--overlap must be non-negative and smaller than half --tile-size")

    core_size = tile_size - 2 * overlap
    index = 0
    row = 0
    for core_y0 in range(0, height, core_size):
        core_y1 = min(core_y0 + core_size, height)
        column = 0
        for core_x0 in range(0, width, core_size):
            core_x1 = min(core_x0 + core_size, width)
            yield Tile(
                index=index,
                row=row,
                column=column,
                core_x0=core_x0,
                core_y0=core_y0,
                core_x1=core_x1,
                core_y1=core_y1,
                read_x0=max(0, core_x0 - overlap),
                read_y0=max(0, core_y0 - overlap),
                read_x1=min(width, core_x1 + overlap),
                read_y1=min(height, core_y1 + overlap),
            )
            index += 1
            column += 1
        row += 1


def _fallback_points(labels: np.ndarray) -> np.ndarray:
    from scipy import ndimage

    count = int(labels.max())
    if count == 0:
        return np.empty((0, 2), dtype=np.float64)
    points = ndimage.center_of_mass(
        np.ones(labels.shape, dtype=np.uint8),
        labels,
        np.arange(1, count + 1),
    )
    return np.asarray(points, dtype=np.float64)


def keep_core_owned_instances(
    labels: np.ndarray,
    points_yx: np.ndarray | None,
    tile: Tile,
) -> tuple[np.ndarray, int]:
    """Keep complete instances whose centroid belongs to this tile's core."""
    labels = np.asarray(labels)
    label_count = int(labels.max())
    if label_count == 0:
        return np.zeros(labels.shape, dtype=np.uint32), 0

    points = np.asarray(points_yx) if points_yx is not None else np.empty((0, 2))
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < label_count:
        points = _fallback_points(labels)

    global_y = points[:label_count, 0] + tile.read_y0
    global_x = points[:label_count, 1] + tile.read_x0
    owned = (
        (global_x >= tile.core_x0)
        & (global_x < tile.core_x1)
        & (global_y >= tile.core_y0)
        & (global_y < tile.core_y1)
    )
    lookup = np.zeros(label_count + 1, dtype=np.uint32)
    owned_ids = np.flatnonzero(owned) + 1
    lookup[owned_ids] = owned_ids.astype(np.uint32)
    return lookup[labels], int(len(owned_ids))


def stitch_tile(
    destination: np.ndarray,
    local_labels: np.ndarray,
    next_label: int,
    merge_overlap: float,
) -> tuple[int, int, int]:
    """Paste a tile, merging near-identical instances in overlapping halos."""
    local = np.asarray(local_labels, dtype=np.uint32)
    dest = np.asarray(destination)
    if local.shape != dest.shape:
        raise ValueError(f"Tile shape mismatch: labels={local.shape}, destination={dest.shape}")
    if not 0.0 <= merge_overlap <= 1.0:
        raise ValueError("--merge-overlap must be between 0 and 1")

    local_ids = np.unique(local)
    local_ids = local_ids[local_ids != 0]
    if local_ids.size == 0:
        return next_label, 0, 0

    local_max = int(local_ids[-1])
    lookup = np.zeros(local_max + 1, dtype=np.uint32)
    local_areas = np.bincount(local.ravel(), minlength=local_max + 1)
    merged_local: set[int] = set()

    conflict = (local != 0) & (dest != 0)
    if np.any(conflict):
        encoded = (local[conflict].astype(np.uint64) << np.uint64(32)) | dest[conflict].astype(np.uint64)
        pairs, counts = np.unique(encoded, return_counts=True)
        best: dict[int, tuple[int, int]] = {}
        for encoded_pair, count in zip(pairs.tolist(), counts.tolist()):
            local_id = int(encoded_pair >> 32)
            global_id = int(encoded_pair & 0xFFFFFFFF)
            previous = best.get(local_id)
            if previous is None or count > previous[1]:
                best[local_id] = (global_id, int(count))
        for local_id, (global_id, overlap_pixels) in best.items():
            if overlap_pixels / max(1, int(local_areas[local_id])) >= merge_overlap:
                lookup[local_id] = global_id
                merged_local.add(local_id)

    fresh_ids = [int(value) for value in local_ids if int(value) not in merged_local]
    if next_label + len(fresh_ids) > np.iinfo(np.uint32).max:
        raise OverflowError("More instances than uint32 labels can represent")
    for local_id in fresh_ids:
        next_label += 1
        lookup[local_id] = next_label

    mapped = lookup[local]
    empty = (mapped != 0) & (dest == 0)
    dest[empty] = mapped[empty]
    return next_label, len(fresh_ids), len(merged_local)


def parse_rectangle(value: str, width: int, height: int) -> tuple[int, int, int, int]:
    try:
        x0, y0, x1, y1 = (int(part.strip()) for part in value.split(","))
    except (TypeError, ValueError):
        raise ValueError(f"Invalid rectangle {value!r}; expected x0,y0,x1,y1") from None
    if not (0 <= x0 < x1 <= width and 0 <= y0 < y1 <= height):
        raise ValueError(
            f"Rectangle {value!r} is outside the {width}x{height} image or has zero area"
        )
    return x0, y0, x1, y1


def exclude_and_compact_labels(
    labels: np.ndarray,
    rectangles: list[tuple[int, int, int, int]],
    chunk_rows: int = 1024,
) -> np.ndarray:
    """Remove instances touching excluded rectangles and compact remaining IDs."""
    if not rectangles:
        return np.empty(0, dtype=np.uint32)
    excluded = []
    for x0, y0, x1, y1 in rectangles:
        excluded.append(np.unique(labels[y0:y1, x0:x1]))
    excluded_ids = np.unique(np.concatenate(excluded))
    excluded_ids = excluded_ids[excluded_ids != 0].astype(np.uint32, copy=False)
    if excluded_ids.size == 0:
        return excluded_ids

    max_label = int(labels.max())
    keep = np.ones(max_label + 1, dtype=bool)
    keep[0] = False
    keep[excluded_ids] = False
    lookup = np.zeros(max_label + 1, dtype=np.uint32)
    compacted = np.cumsum(keep, dtype=np.uint32)
    lookup[keep] = compacted[keep]
    for y0 in range(0, labels.shape[0], chunk_rows):
        labels[y0 : y0 + chunk_rows] = lookup[labels[y0 : y0 + chunk_rows]]
    return excluded_ids


def _prepare_tile(tile_image: Image.Image, input_channels: int, channel: str) -> np.ndarray:
    if input_channels > 1:
        return np.asarray(tile_image.convert("RGB"), dtype=np.uint8)
    if channel == "gray":
        selected = tile_image.convert("L")
    else:
        index = {"red": 0, "green": 1, "blue": 2}[channel]
        selected = tile_image.convert("RGB").getchannel(index)
    return np.asarray(selected, dtype=np.uint8)


def _load_model(model_name: str, model_dir: str | None):
    from stardist.models import StarDist2D

    if model_dir:
        model_path = Path(model_dir).expanduser().resolve()
        for required in ("config.json", "thresholds.json", "weights_best.h5"):
            if not (model_path / required).exists():
                raise FileNotFoundError(f"Local StarDist model is missing {required}: {model_path / required}")
        log(f"Loading StarDist model from local directory: {model_path}")
        return StarDist2D(None, name=model_path.name, basedir=str(model_path.parent)), str(model_path)
    log(f"Loading StarDist model: {model_name}")
    return StarDist2D.from_pretrained(model_name), None


def normalize_tile(image: np.ndarray, pmin: float, pmax: float) -> tuple[np.ndarray, str]:
    """Keep ordinary percentile normalization; handle degenerate channels safely.

    A mostly black fluorescence tile can have both percentiles equal to zero
    despite containing real nuclei. Dividing by csbdeep's epsilon (1e-20)
    sends those pixels to ~1e22 and can abort native polygon NMS. Fall back
    to the actual per-channel range, without dropping the sparse signal.
    """
    from csbdeep.utils import normalize_mi_ma

    if not (0 <= pmin < pmax <= 100):
        raise ValueError("Normalization percentiles must satisfy 0 <= low < high <= 100.")
    if image.size == 0 or not np.isfinite(image).all():
        raise ValueError("StarDist tile must contain finite, nonempty image data.")
    if not np.any(image):
        return np.zeros(image.shape, dtype=np.float32), "blank"
    low = np.percentile(image, pmin, axis=(0, 1), keepdims=True)
    high = np.percentile(image, pmax, axis=(0, 1), keepdims=True)
    degenerate = high <= low
    mode = "percentile"
    if np.any(degenerate):
        mode = "range_fallback"
        low = np.where(degenerate, np.min(image, axis=(0, 1), keepdims=True), low)
        high = np.where(degenerate, np.max(image, axis=(0, 1), keepdims=True), high)
        # Constant channels normalize to zero, with a nonzero denominator.
        high = np.where(high <= low, low + 1.0, high)
    normalized = normalize_mi_ma(image, low, high)
    if not np.isfinite(normalized).all():
        raise ValueError("Nonfinite normalized tile; refusing to pass it to native StarDist NMS.")
    return normalized, mode


def predict_tile(model, image: np.ndarray, args) -> tuple[np.ndarray, dict, str]:
    normalized, mode = normalize_tile(image, float(args.normalization_low), float(args.normalization_high))
    if mode == "blank":
        return np.zeros(image.shape[:2], dtype=np.uint32), {"points": np.empty((0, 2))}, mode
    kwargs: dict[str, Any] = {
        "prob_thresh": float(args.prob_thresh),
        "n_tiles": parse_tiles(args.n_tiles, normalized.ndim),
    }
    if args.nms_thresh is not None:
        kwargs["nms_thresh"] = float(args.nms_thresh)
    if float(args.scale) != 1.0:
        kwargs["scale"] = float(args.scale)
    labels, details = model.predict_instances(normalized, **kwargs)
    return labels, details, mode


def _legacy_unsafe_tiles(image_path, tiles, cache_dir, input_channels, channel, pmin, pmax):
    """Only reuse legacy tiles whose normalization is unchanged by the fix."""
    unsafe = set()
    with ImageRegionSource(image_path) as source:
        for tile in tiles:
            if not _tile_cache_path(cache_dir, tile).is_file():
                continue
            pixels = _prepare_tile(source.crop(tile.read_box), input_channels, channel)
            _, mode = normalize_tile(pixels, pmin, pmax)
            if mode != "percentile":
                unsafe.add(tile.index)
    return unsafe


def _tile_cache_path(cache_dir: Path, tile: Tile) -> Path:
    return cache_dir / f"tile_r{tile.row:03d}_c{tile.column:03d}.owned_labels.tif"


def _fingerprint(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _write_preview(
    image_path: Path,
    labels_path: Path,
    output_path: Path,
    max_dim: int,
) -> None:
    if max_dim <= 0:
        return
    with ImageRegionSource(image_path) as source:
        width, height = source.size
        scale = min(1.0, float(max_dim) / max(width, height))
        preview_size = (max(1, round(width * scale)), max(1, round(height * scale)))
        preview = source.preview_rgb(max_dim)

    labels = tifffile.memmap(labels_path, mode="r")
    y_index = np.minimum(
        (np.arange(preview_size[1], dtype=np.float64) / scale).astype(np.int64),
        height - 1,
    )
    x_index = np.minimum(
        (np.arange(preview_size[0], dtype=np.float64) / scale).astype(np.int64),
        width - 1,
    )
    sampled = np.asarray(labels[np.ix_(y_index, x_index)])
    from skimage.segmentation import find_boundaries

    preview_array = np.asarray(preview).copy()
    preview_array[find_boundaries(sampled, mode="outer")] = np.array([255, 40, 40], dtype=np.uint8)
    Image.fromarray(preview_array).save(output_path, quality=92)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run full-resolution StarDist inference using overlapping disk-cached tiles."
    )
    parser.add_argument("--image", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--model", default=None)
    parser.add_argument("--model-dir")
    parser.add_argument("--prob-thresh", type=float, default=0.30)
    parser.add_argument("--nms-thresh", type=float, default=None)
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument("--n-tiles", default="2,2", help="StarDist compute tiling within each disk tile.")
    parser.add_argument("--tile-size", type=int, default=4096)
    parser.add_argument("--overlap", type=int, default=256)
    parser.add_argument("--merge-overlap", type=float, default=0.50)
    parser.add_argument("--channel", choices=("gray", "red", "green", "blue"), default="gray")
    parser.add_argument("--normalization-low", type=float, default=1.0)
    parser.add_argument("--normalization-high", type=float, default=99.8)
    parser.add_argument("--preview-max-dim", type=int, default=3000)
    parser.add_argument("--cache-dir", help="Tile cache directory. Defaults to <outdir>/tile_cache.")
    parser.add_argument(
        "--exclude-rectangle",
        action="append",
        default=[],
        metavar="X0,Y0,X1,Y1",
        help="Remove instances touching this image rectangle after reconstruction; may be repeated.",
    )
    parser.add_argument("--no-resume", action="store_true", help="Reject existing cached tiles.")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    if not (0 <= args.normalization_low < args.normalization_high <= 100):
        raise ValueError("Normalization percentiles must satisfy 0 <= low < high <= 100.")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

    image_path = Path(args.image).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    cache_dir = Path(args.cache_dir).expanduser().resolve() if args.cache_dir else outdir / "tile_cache"
    outdir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    if not image_path.is_file():
        raise FileNotFoundError(f"Input image not found: {image_path}")

    with ImageRegionSource(image_path) as source:
        width, height = source.size
    tiles = list(iter_tiles(width, height, int(args.tile_size), int(args.overlap)))
    model_name = args.model or _default_model_name(image_path)
    model_dir = args.model_dir or _default_model_dir(model_name)
    image_stat = image_path.stat()
    run_parameters = {
        "image": str(image_path),
        "image_size": [width, height],
        "image_bytes": image_stat.st_size,
        "image_mtime_ns": image_stat.st_mtime_ns,
        "model": model_name,
        "model_dir": str(Path(model_dir).expanduser().resolve()) if model_dir else None,
        "prob_thresh": float(args.prob_thresh),
        "nms_thresh": args.nms_thresh,
        "scale": float(args.scale),
        "n_tiles": args.n_tiles,
        "tile_size": int(args.tile_size),
        "overlap": int(args.overlap),
        "channel": args.channel,
        "normalization_low": float(args.normalization_low),
        "normalization_high": float(args.normalization_high),
    }
    manifest_path = cache_dir / "manifest.json"
    fingerprint = _fingerprint(run_parameters)
    legacy_cache = False
    if manifest_path.exists():
        old_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if old_manifest.get("fingerprint") != fingerprint:
            raise RuntimeError(
                f"Tile cache parameters do not match this run: {manifest_path}. "
                "Use a different --cache-dir or remove the old cache."
            )
        if args.no_resume and any(_tile_cache_path(cache_dir, tile).exists() for tile in tiles):
            raise RuntimeError(f"Cached tiles exist but --no-resume was requested: {cache_dir}")
        policy = old_manifest.get("normalization_policy")
        if policy not in (None, NORMALIZATION_POLICY):
            raise RuntimeError(f"Unsupported normalization policy in tile cache: {policy}")
        legacy_cache = policy is None
    else:
        if any(_tile_cache_path(cache_dir, tile).exists() for tile in tiles):
            raise RuntimeError(
                f"Cached tiles exist without a manifest: {cache_dir}. "
                "Cannot verify their image/parameters; use a different --cache-dir."
            )
        manifest_path.write_text(
            json.dumps(
                {
                    "fingerprint": fingerprint,
                    "parameters": run_parameters,
                    "tiles": [asdict(tile) for tile in tiles],
                    "normalization_policy": NORMALIZATION_POLICY,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )

    missing_tiles = [tile for tile in tiles if not _tile_cache_path(cache_dir, tile).is_file()]
    model = None
    resolved_model_dir = str(Path(model_dir).expanduser().resolve()) if model_dir else None
    input_channels = None
    inference_start = time.time()
    if missing_tiles or legacy_cache:
        model, resolved_model_dir = _load_model(model_name, model_dir)
        input_channels = int(getattr(model.config, "n_channel_in", 1) or 1)
    if legacy_cache:
        unsafe = _legacy_unsafe_tiles(
            image_path, tiles, cache_dir, input_channels, args.channel,
            float(args.normalization_low), float(args.normalization_high),
        )
        missing_tiles = [tile for tile in tiles
                         if tile.index in unsafe or not _tile_cache_path(cache_dir, tile).is_file()]
        log(f"Legacy cache checked: reuse={len(tiles) - len(missing_tiles)}, "
            f"recompute_degenerate={len(unsafe)}; policy={NORMALIZATION_POLICY}")
    log(
        f"image_size={width}x{height}, disk_tiles={len(tiles)}, missing={len(missing_tiles)}, "
        f"tile_size={args.tile_size}, overlap={args.overlap}, model={model_name}, channel={args.channel}"
    )
    normalization_events = []
    if missing_tiles:
        with ImageRegionSource(image_path) as source:
            log(f"image_region_backend={source.backend}")
            for completed, tile in enumerate(missing_tiles, start=1):
                started = time.time()
                tile_image = source.crop(tile.read_box)
                image_array = _prepare_tile(tile_image, input_channels, args.channel)
                labels, details, normalization_mode = predict_tile(model, image_array, args)
                if normalization_mode != "percentile":
                    normalization_events.append({"tile": tile.index + 1, "mode": normalization_mode})
                if labels.shape != tile.read_shape:
                    raise RuntimeError(
                        f"StarDist returned {labels.shape} for tile {tile.index}; expected {tile.read_shape}"
                    )
                owned, owned_count = keep_core_owned_instances(labels, details.get("points"), tile)
                cache_path = _tile_cache_path(cache_dir, tile)
                temporary = cache_path.with_suffix(".tmp.tif")
                tifffile.imwrite(temporary, owned, compression="zlib", predictor=True)
                os.replace(temporary, cache_path)
                log(
                    f"tile {tile.index + 1}/{len(tiles)} (missing {completed}/{len(missing_tiles)}): "
                    f"read={tile.read_box}, labels={int(labels.max())}, owned={owned_count}, "
                    f"normalization={normalization_mode}, "
                    f"elapsed={time.time() - started:.1f}s"
                )
    del model
    if legacy_cache:
        # Upgrade only once every unsafe/missing tile has been replaced.
        # If inference aborts, the next run still validates all legacy tiles.
        old_manifest["normalization_policy"] = NORMALIZATION_POLICY
        temporary_manifest = manifest_path.with_suffix(".tmp.json")
        temporary_manifest.write_text(json.dumps(old_manifest, indent=2, sort_keys=True), encoding="utf-8")
        os.replace(temporary_manifest, manifest_path)

    labels_path = outdir / f"{args.sample}.stardist_labels.tif"
    temporary_labels = outdir / f".{args.sample}.{os.getpid()}.stardist_labels.tmp.tif"
    log(f"Rebuilding full-resolution labels: {labels_path}")
    global_labels = tifffile.memmap(
        temporary_labels,
        shape=(height, width),
        dtype=np.uint32,
        bigtiff=True,
        photometric="minisblack",
    )
    next_label = 0
    total_fresh = 0
    total_merged = 0
    for tile in tiles:
        local = tifffile.imread(_tile_cache_path(cache_dir, tile))
        destination = global_labels[tile.read_y0 : tile.read_y1, tile.read_x0 : tile.read_x1]
        next_label, fresh, merged = stitch_tile(
            destination,
            local,
            next_label,
            float(args.merge_overlap),
        )
        total_fresh += fresh
        total_merged += merged
    excluded_rectangles = [
        parse_rectangle(value, width, height) for value in args.exclude_rectangle
    ]
    excluded_ids = exclude_and_compact_labels(global_labels, excluded_rectangles)
    final_label_count = next_label - int(len(excluded_ids))
    global_labels.flush()
    del global_labels
    os.replace(temporary_labels, labels_path)

    bbox_path = outdir / f"{args.sample}.stardist_bbox.csv"
    bbox_path.write_text(f"0,0,{width},{height}\n", encoding="utf-8")
    preview_path = outdir / f"{args.sample}.stardist_boundary_overlay.preview.jpg"
    _write_preview(image_path, labels_path, preview_path, int(args.preview_max_dim))
    summary = {
        "sample": args.sample,
        "image": str(image_path),
        "model": model_name,
        "model_dir": resolved_model_dir,
        "input_channels": input_channels,
        "channel": args.channel,
        "prob_thresh": float(args.prob_thresh),
        "nms_thresh": args.nms_thresh,
        "scale": float(args.scale),
        "n_tiles": args.n_tiles,
        "tile_size": int(args.tile_size),
        "overlap": int(args.overlap),
        "merge_overlap": float(args.merge_overlap),
        "disk_tiles": len(tiles),
        "original_size": [width, height],
        "processed_size": [width, height],
        "resize_scale": 1.0,
        "labels": final_label_count,
        "fresh_instances": total_fresh,
        "merged_instances": total_merged,
        "excluded_instances": int(len(excluded_ids)),
        "excluded_rectangles": [list(rectangle) for rectangle in excluded_rectangles],
        "labels_path": str(labels_path),
        "latest_labels": str(labels_path),
        "bbox": str(bbox_path),
        "boundary_overlay_preview": str(preview_path),
        "tile_cache": str(cache_dir),
        "normalization_policy": NORMALIZATION_POLICY,
        "normalization_events_this_run": normalization_events,
        "reused_tiles": len(tiles) - len(missing_tiles),
        "inference_elapsed_min": (time.time() - inference_start) / 60.0,
    }
    summary_path = outdir / f"{args.sample}.stardist_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    log(
        f"labels={final_label_count:,}, merged={total_merged:,}, "
        f"excluded={len(excluded_ids):,}, labels_path={labels_path}"
    )
    log(f"summary={summary_path}")


if __name__ == "__main__":
    main()
