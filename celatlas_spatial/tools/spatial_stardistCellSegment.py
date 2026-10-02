"""Step 06.segment/02: StarDist nucleus segmentation and cell-level assignment."""

from __future__ import annotations

import json
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path
import importlib.util

from celatlas_spatial.rna.mkref import Mkref_rna


DEFAULT_ENV_NAME = "celatlas18"


def _workspace_path(*parts: str) -> Path:
    workspace = os.environ.get("CELATLAS_WORKSPACE")
    root = Path(workspace).expanduser() if workspace else Path.home() / "celatlas_spatial"
    return root.joinpath(*parts)


def _is_true(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _infer_sample(args) -> str:
    if args.sample:
        return args.sample
    if args.sampledir:
        return Path(args.sampledir).resolve().name
    raise ValueError("Cannot infer sample. Provide --sample or --sampledir.")


def _find_registered_image(bin_segment_dir: Path, sample: str) -> Path:
    candidates = [
        bin_segment_dir / "images" / f"{sample}_regist.tif",
        bin_segment_dir / "images" / f"{sample}_regist.tiff",
        bin_segment_dir / "images" / f"{sample}_regist.png",
        bin_segment_dir / "images" / f"{sample}_regist.jpg",
        bin_segment_dir / "images" / f"{sample}_regist.jpeg",
        bin_segment_dir / "images" / "tissue_hires_image.png",
        bin_segment_dir / "images" / "tissue_hires_image.jpg",
        bin_segment_dir / "images" / "tissue_hires_image.jpeg",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find registered tissue image. Provide --image. Tried: "
        + ", ".join(str(path) for path in candidates)
    )


def _find_barcode_mask(bin_segment_dir: Path, sample: str, *, include_tissue: bool = False) -> Path | None:
    mask_dir = bin_segment_dir / "cavity_filter" / "mask"
    candidates = [
        mask_dir / f"{sample}_tissue_mask.png",
        mask_dir / f"{sample}_regist_tissue_mask.png",
    ]
    if include_tissue:
        # This opt-in is used by the DAPI preset: older HE cavity masks
        # are not applicable to a newly registered fluorescence image.
        candidates = [bin_segment_dir / "images" / f"{sample}_tissue_cut.tif"]
    for path in candidates:
        if path.exists():
            return path.resolve()
    return None


FLUORESCENCE_NAME_HINTS = ("dapi", "fluores", "fluoro", "immunofluores", "multiplex", "mif")


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


def _env_python(path_value: str | None) -> str | None:
    if not path_value:
        return None
    python_path = Path(path_value).expanduser() / "bin" / "python"
    if python_path.exists():
        return str(python_path)
    return None


def _default_stardist_python() -> str:
    for env_name in ("STARDIST_PYTHON", "CELATLAS_STARDIST_PYTHON"):
        env_value = os.environ.get(env_name)
        if env_value:
            return env_value

    env_value = _env_python(os.environ.get("CELATLAS_ENV_PREFIX"))
    if env_value:
        return env_value

    if importlib.util.find_spec("stardist") is not None:
        return sys.executable

    env_value = _env_python(os.environ.get("CONDA_PREFIX"))
    if env_value:
        return env_value

    conda_root = os.environ.get("CELATLAS_CONDA_ROOT") or os.environ.get("MAMBA_ROOT_PREFIX")
    env_name = os.environ.get("CELATLAS_ENV_NAME", DEFAULT_ENV_NAME)
    if conda_root:
        env_value = _env_python(str(Path(conda_root).expanduser() / "envs" / env_name))
        if env_value:
            return env_value

    return sys.executable


def _resolve_paths(args) -> dict[str, Path | str | None]:
    sample = _infer_sample(args)
    sampledir = Path(args.sampledir).resolve() if args.sampledir else None
    if sampledir is None:
        sampledir = Path(os.environ.get("CELATLAS_RESULTS_ROOT", _workspace_path("results"))) / sample

    if args.bin_segment_dir:
        bin_segment_dir = Path(args.bin_segment_dir).resolve()
    else:
        bin_segment_dir = sampledir / "06.segment" / "01.binsegment"
    runtime_input_dir = (
        Path(args.runtime_input_dir).resolve()
        if args.runtime_input_dir
        else sampledir / "06.segment" / "mask"
    )
    outdir = Path(args.outdir).resolve() if args.outdir else sampledir / "06.segment" / "02.cellsegment"
    image_path = Path(args.image).resolve() if args.image else _find_registered_image(bin_segment_dir, sample)
    barcode_positions = (
        Path(args.barcode_positions).resolve()
        if args.barcode_positions
        else runtime_input_dir / f"{sample}_FilterBarcodes.csv"
    )
    count_detail = (
        Path(args.count_detail).resolve()
        if args.count_detail
        else sampledir / "05.count" / f"{sample}_count_detail.txt"
    )
    tissue_bbox = (
        Path(args.tissue_bbox).resolve()
        if args.tissue_bbox
        else runtime_input_dir / f"{sample}_tissue_bbox.csv"
    )
    labels = Path(args.labels).resolve() if args.labels else None
    barcode_mask = Path(args.barcode_mask).resolve() if args.barcode_mask else None
    if barcode_mask is None:
        barcode_mask = _find_barcode_mask(
            bin_segment_dir, sample, include_tissue=getattr(args, "require_registered_dapi", False)
        )
    gtf = None
    if args.genomeDir:
        gtf = Path(Mkref_rna.parse_genomeDir(args.genomeDir)["gtf"]).resolve()

    return {
        "sample": sample,
        "sampledir": sampledir,
        "bin_segment_dir": bin_segment_dir,
        "runtime_input_dir": runtime_input_dir,
        "outdir": outdir,
        "image": image_path,
        "barcode_positions": barcode_positions,
        "count_detail": count_detail,
        "tissue_bbox": tissue_bbox,
        "labels": labels,
        "barcode_mask": barcode_mask,
        "gtf": gtf,
    }


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def _validate_registered_dapi(paths: dict) -> None:
    """Reject legacy white-background images and mismatched tissue canvases."""
    import tifffile
    from PIL import Image

    image = Path(paths["image"])
    expected = Path(paths["bin_segment_dir"]) / "images" / f"{paths['sample']}_regist.tif"
    if image.resolve() != expected.resolve():
        raise ValueError("DAPI preset requires the full-resolution binSegment registered TIFF.")
    _require_file(image, "Full-resolution registered DAPI image")
    provenance = image.with_name(image.name + ".registration.json")
    hint = "Rerun binSegment with --cell-segmentation-preset dapi; do not use --skip-binsegment with legacy white-background images."
    if not provenance.is_file():
        raise ValueError(f"Registered DAPI provenance is missing. {hint}")
    metadata = json.loads(provenance.read_text(encoding="utf-8"))
    stat = image.stat()
    with tifffile.TiffFile(image) as tif:
        shape = tuple(tif.series[0].shape)
        valid_layout = len(tif.pages) == 1 and tif.series[0].axes in {"YX", "YXS"}
    if (not valid_layout or metadata.get("method") != "HE"
            or metadata.get("fluorescence_background") is not True
            or metadata.get("shape") != list(shape)
            or metadata.get("image_bytes") != stat.st_size
            or metadata.get("image_mtime_ns") != stat.st_mtime_ns):
        raise ValueError(f"Registered image is not a verified black-background DAPI canvas. {hint}")
    if paths["barcode_mask"] is None:
        raise ValueError(f"Final tissue mask is missing. {hint}")
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(paths["barcode_mask"]) as mask:
        if mask.size != (shape[1], shape[0]):
            raise ValueError("DAPI barcode mask must match the full-resolution registered image canvas.")


def _run_command(cmd: list[str], env: dict[str, str] | None = None) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _resolve_assignment_bbox(args, paths: dict[str, Path | str | None]) -> tuple[Path | None, str]:
    """Resolve the full-resolution coordinate canvas used by StarDist labels."""
    if args.tissue_bbox:
        bbox = Path(args.tissue_bbox).resolve()
        _require_file(bbox, "Explicit tissue bbox")
        return bbox, "explicit"

    sample = str(paths["sample"])
    outdir = Path(paths["outdir"])
    inference_bbox = outdir / "stardist_inference" / f"{sample}.stardist_bbox.csv"
    if inference_bbox.exists():
        return inference_bbox, "stardist_inference"

    image_path = Path(paths["image"])
    if image_path.exists():
        from PIL import Image

        Image.MAX_IMAGE_PIXELS = None
        with Image.open(image_path) as image:
            width, height = image.size
        inference_bbox.parent.mkdir(parents=True, exist_ok=True)
        inference_bbox.write_text(f"0,0,{width},{height}\n", encoding="utf-8")
        return inference_bbox, "registered_image"

    runtime_bbox = Path(paths["tissue_bbox"])
    if runtime_bbox.exists():
        return runtime_bbox, "runtime_input"
    return None, "label_shape"


def _run_inference(args, paths: dict[str, Path | str | None]) -> Path:
    sample = str(paths["sample"])
    outdir = Path(paths["outdir"])
    inference_dir = outdir / "stardist_inference"
    inference_dir.mkdir(parents=True, exist_ok=True)

    stardist_python = args.stardist_python or _default_stardist_python()
    python_path = Path(stardist_python)
    if not python_path.exists() or not os.access(python_path, os.X_OK):
        raise FileNotFoundError(f"StarDist Python is not executable: {stardist_python}")

    image_path = Path(paths["image"])
    _require_file(image_path, "StarDist input image")
    selected_model = args.model or _default_model_name(image_path)

    inference_script = Path(__file__).resolve().with_name(
        "stardist_tiled_inference.py" if args.tiled_inference else "stardist_inference.py"
    )
    cmd = [
        str(python_path),
        str(inference_script),
        "--image",
        str(image_path),
        "--outdir",
        str(inference_dir),
        "--sample",
        sample,
        "--model",
        selected_model,
        "--prob-thresh",
        str(args.prob_thresh),
        "--scale",
        str(args.scale),
        "--n-tiles",
        args.n_tiles,
    ]
    if args.tiled_inference:
        cmd.extend(
            [
                "--tile-size",
                str(args.tile_size),
                "--overlap",
                str(args.tile_overlap),
                "--merge-overlap",
                str(args.tile_merge_overlap),
                "--channel",
                args.fluorescence_channel,
            ]
        )
        if args.tile_cache_dir:
            cmd.extend(["--cache-dir", args.tile_cache_dir])
        for rectangle in args.exclude_rectangle:
            cmd.extend(["--exclude-rectangle", rectangle])
    else:
        cmd.extend(["--max-dim", str(args.max_dim)])
    model_dir = args.model_dir or _default_model_dir(selected_model)
    if model_dir:
        cmd.extend(["--model-dir", model_dir])
    if args.nms_thresh is not None:
        cmd.extend(["--nms-thresh", str(args.nms_thresh)])
    if args.preview_max_dim is not None:
        cmd.extend(["--preview-max-dim", str(args.preview_max_dim)])

    if getattr(args, "require_registered_dapi", False):
        if model_dir:
            model_config = json.loads((Path(model_dir) / "config.json").read_text(encoding="utf-8"))
            if int(model_config.get("n_channel_in", 1)) != 1:
                raise ValueError("DAPI requires a single-channel fluorescence model; check STARDIST_MODEL_DIR / --model-dir.")
        if not args.tile_cache_dir:
            # Re-registration changes the canvas/input signature. Keep earlier
            # tiles intact, but never reuse them for a different image/run.
            stat = image_path.stat()
            signature = {"command": cmd, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
            if model_dir:
                signature["model_files"] = [
                    (p.name, p.stat().st_size, p.stat().st_mtime_ns)
                    for p in sorted(Path(model_dir).iterdir()) if p.is_file()
                ]
            key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:20]
            args.tile_cache_dir = str(inference_dir / "tile_cache" / key)
            cmd.extend(["--cache-dir", args.tile_cache_dir])

    env = os.environ.copy()
    mplconfigdir = args.mplconfigdir or str(outdir / ".mplconfig")
    Path(mplconfigdir).mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = mplconfigdir
    numba_cache_dir = outdir / ".numba_cache"
    numba_cache_dir.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_cache_dir)
    env.setdefault("PYTHONNOUSERSITE", "1")
    stardist_prefix = python_path.resolve().parent.parent
    stardist_lib = stardist_prefix / "lib"
    if stardist_lib.exists():
        current_ld = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            str(stardist_lib) if not current_ld else f"{stardist_lib}:{current_ld}"
        )

    _run_command(cmd, env=env)
    labels = inference_dir / f"{sample}.stardist_labels.tif"
    _require_file(labels, "StarDist labels")
    return labels


def _run_assignment(args, paths: dict[str, Path | str | None], labels: Path) -> dict[str, str | None]:
    sample = str(paths["sample"])
    outdir = Path(paths["outdir"])
    barcode_positions = Path(paths["barcode_positions"])
    count_detail = Path(paths["count_detail"])
    image_path = Path(paths["image"])
    barcode_mask = paths["barcode_mask"]
    gtf = paths.get("gtf")

    _require_file(labels, "StarDist labels")
    _require_file(barcode_positions, "Barcode positions")
    if not _is_true(args.skip_counts):
        _require_file(count_detail, "Count detail")

    assignment_python = args.assignment_python or sys.executable
    python_path = Path(assignment_python)
    if not python_path.exists() or not os.access(python_path, os.X_OK):
        raise FileNotFoundError(f"Assignment Python is not executable: {assignment_python}")

    assignment_script = Path(__file__).resolve().with_name("stardist_assign.py")
    cmd = [
        str(python_path),
        str(assignment_script),
        "--labels",
        str(labels),
        "--barcode-positions",
        str(barcode_positions),
        "--sample",
        sample,
        "--expand-pixels",
        str(args.expand_pixels),
        "--min-umi",
        str(args.min_umi),
        "--min-genes",
        str(args.min_genes),
        "--outdir",
        str(outdir),
        "--rounding",
        args.rounding,
        "--chunksize",
        str(args.chunksize),
        "--visualization-preview-dim",
        str(args.visualization_preview_dim),
    ]
    if not _is_true(args.skip_counts):
        cmd.extend(["--count-detail", str(count_detail)])
    assignment_bbox, bbox_source = _resolve_assignment_bbox(args, paths)
    if assignment_bbox is not None:
        cmd.extend(["--tissue-bbox", str(assignment_bbox)])
    print(
        "Assignment coordinate bbox: "
        f"{assignment_bbox if assignment_bbox is not None else 'label shape'} "
        f"(source={bbox_source})"
    )
    if image_path.exists():
        cmd.extend(["--he-image", str(image_path)])
    if gtf:
        gtf_path = Path(gtf)
        _require_file(gtf_path, "GTF")
        cmd.extend(["--gtf", str(gtf_path)])
    if args.max_umi is not None:
        cmd.extend(["--max-umi", str(args.max_umi)])
    for allowlist in args.barcode_allowlist or []:
        _require_file(Path(allowlist), "Barcode allowlist")
        cmd.extend(["--barcode-allowlist", allowlist])
    if barcode_mask:
        barcode_mask_path = Path(barcode_mask)
        _require_file(barcode_mask_path, "Barcode mask")
        cmd.extend(
            [
                "--barcode-mask",
                str(barcode_mask_path),
                "--barcode-mask-keep-value",
                args.barcode_mask_keep_value,
                "--barcode-mask-threshold",
                str(args.barcode_mask_threshold),
            ]
        )
    if args.barcode_in_tissue_only:
        cmd.append("--barcode-in-tissue-only")
    if args.no_label_output:
        cmd.append("--no-label-output")

    env = os.environ.copy()
    mplconfigdir = args.mplconfigdir or str(outdir / ".mplconfig")
    Path(mplconfigdir).mkdir(parents=True, exist_ok=True)
    env["MPLCONFIGDIR"] = mplconfigdir
    numba_cache_dir = outdir / ".numba_cache"
    numba_cache_dir.mkdir(parents=True, exist_ok=True)
    env["NUMBA_CACHE_DIR"] = str(numba_cache_dir)
    env.setdefault("PYTHONNOUSERSITE", "1")
    assignment_prefix = python_path.resolve().parent.parent
    assignment_lib = assignment_prefix / "lib"
    if assignment_lib.exists():
        current_ld = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = (
            str(assignment_lib) if not current_ld else f"{assignment_lib}:{current_ld}"
        )
    _run_command(cmd, env=env)
    return {
        "path": str(assignment_bbox) if assignment_bbox is not None else None,
        "source": bbox_source,
    }


def run_stardist_cell_segment(args) -> dict:
    paths = _resolve_paths(args)
    if getattr(args, "require_registered_dapi", False):
        _validate_registered_dapi(paths)
    sample = str(paths["sample"])
    outdir = Path(paths["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Sample: {sample}")
    print(f"Sample dir: {paths['sampledir']}")
    print(f"Output dir: {outdir}")
    print(f"Registered tissue image: {paths['image']}")
    print(f"Barcode positions: {paths['barcode_positions']}")
    print(f"Count detail: {paths['count_detail']}")
    print(f"Barcode mask: {paths['barcode_mask'] or 'None'}")

    labels = Path(paths["labels"]) if paths["labels"] else None
    if labels is None:
        if args.skip_inference:
            raise ValueError("--skip-inference requires --labels")
        labels = _run_inference(args, paths)
    else:
        _require_file(labels, "StarDist labels")

    assignment_bbox = _run_assignment(args, paths, labels)
    selected_model = args.model or _default_model_name(Path(paths["image"]))
    resolved_model_dir = args.model_dir or _default_model_dir(selected_model)

    run_summary = {
        "sample": sample,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "sampledir": str(paths["sampledir"]),
        "bin_segment_dir": str(paths["bin_segment_dir"]),
        "runtime_input_dir": str(paths["runtime_input_dir"]),
        "outdir": str(outdir),
        "image": str(paths["image"]),
        "labels": str(labels),
        "barcode_positions": str(paths["barcode_positions"]),
        "count_detail": str(paths["count_detail"]),
        "tissue_bbox": str(paths["tissue_bbox"]),
        "assignment_bbox": assignment_bbox,
        "barcode_mask": str(paths["barcode_mask"]) if paths["barcode_mask"] else None,
        "gtf": str(paths["gtf"]) if paths.get("gtf") else None,
        "skip_counts": bool(args.skip_counts),
        "stardist": {
            "model": selected_model,
            "model_dir": resolved_model_dir,
            "prob_thresh": float(args.prob_thresh),
            "nms_thresh": args.nms_thresh,
            "max_dim": int(args.max_dim),
            "scale": float(args.scale),
            "n_tiles": args.n_tiles,
            "tiled_inference": bool(args.tiled_inference),
            "tile_size": int(args.tile_size),
            "tile_overlap": int(args.tile_overlap),
            "tile_merge_overlap": float(args.tile_merge_overlap),
            "fluorescence_channel": args.fluorescence_channel,
            "tile_cache_dir": args.tile_cache_dir,
        },
        "assignment": {
            "expand_pixels": float(args.expand_pixels),
            "min_umi": int(args.min_umi),
            "min_genes": int(args.min_genes),
            "max_umi": args.max_umi,
        },
    }
    _write_json(outdir / f"{sample}.stardist_cell_segment_run_summary.json", run_summary)
    print(f"StarDist cell segment completed: {outdir}")
    return run_summary


def stardistCellSegment(args):
    run_stardist_cell_segment(args)


def get_opts_stardistCellSegment(parser, sub_program):
    parser.add_argument("--sampledir", help="Sample output directory containing 05.count and 06.segment/01.binsegment.")
    parser.add_argument("--sample", help="Sample name. Defaults to basename of --sampledir.")
    parser.add_argument("--genomeDir", help="Reference genome directory; used to map gene_id to gene_name for cell matrices.")
    parser.add_argument(
        "--runtime-input-dir",
        help="Runtime mask/input directory containing FilterBarcodes and tissue_bbox CSV; defaults to <sampledir>/06.segment/mask.",
    )
    parser.add_argument(
        "--bin-segment-dir",
        help="Path to binSegment output. Defaults to <sampledir>/06.segment/01.binsegment.",
    )
    parser.add_argument("--image", help="Registered tissue image. Defaults to <bin-segment-dir>/images/<sample>_regist.tif.")
    parser.add_argument("--outdir", help="Output directory. Defaults to <sampledir>/06.segment/02.cellsegment.")
    parser.add_argument("--labels", help="Existing StarDist label TIFF. If omitted, inference is run.")
    parser.add_argument("--skip-inference", action="store_true", help="Require --labels and skip StarDist inference.")
    parser.add_argument("--skip-counts", action="store_true", help="Run inference/barcode assignment only; do not aggregate count_detail.")
    parser.add_argument("--require-registered-dapi", action="store_true",
                        help="Require a verified black-background registered TIFF and matching final tissue mask (DAPI runner preset).")

    parser.add_argument("--stardist-python", help="Python executable with stardist/tensorflow installed.")
    parser.add_argument("--assignment-python", help="Python executable for barcode/count assignment. Defaults to current Python.")
    parser.add_argument(
        "--model",
        default=None,
        help=(
            "StarDist pretrained model name. Defaults to 2D_versatile_fluo for fluorescence-like "
            "inputs and 2D_versatile_he otherwise."
        ),
    )
    parser.add_argument("--model-dir", help="Local StarDist model directory; defaults to STARDIST_MODEL_DIR or local cache if present.")
    parser.add_argument("--prob-thresh", type=float, default=0.30, help="StarDist probability threshold.")
    parser.add_argument("--nms-thresh", type=float, default=None, help="Optional StarDist NMS threshold.")
    parser.add_argument("--max-dim", type=int, default=6000, help="Resize image long edge before StarDist; 0 disables.")
    parser.add_argument("--scale", type=float, default=2.0, help="StarDist internal spatial scale.")
    parser.add_argument("--n-tiles", default="4,4", help="Comma-separated StarDist tile grid, e.g. 4,4.")
    parser.add_argument("--preview-max-dim", type=int, default=3000, help="Long edge for StarDist boundary preview.")
    parser.add_argument(
        "--tiled-inference",
        action="store_true",
        help="Run full-resolution disk-cached overlapping-tile inference instead of resizing the whole image.",
    )
    parser.add_argument("--tile-size", type=int, default=4096, help="Input tile width and height for tiled inference.")
    parser.add_argument("--tile-overlap", type=int, default=256, help="Context halo around each tiled-inference core.")
    parser.add_argument(
        "--tile-merge-overlap",
        type=float,
        default=0.50,
        help="Minimum fraction of a tile instance overlapping an existing instance to merge them.",
    )
    parser.add_argument(
        "--fluorescence-channel",
        choices=("gray", "red", "green", "blue"),
        default="gray",
        help="RGB channel used by a one-channel StarDist model.",
    )
    parser.add_argument("--tile-cache-dir", help="Optional resumable tiled-inference cache directory.")
    parser.add_argument(
        "--exclude-rectangle",
        action="append",
        default=[],
        metavar="X0,Y0,X1,Y1",
        help="Remove tile-inference instances touching this annotation/scale-bar rectangle; may be repeated.",
    )
    parser.add_argument("--mplconfigdir", help="Writable matplotlib cache directory for StarDist/visualization.")

    parser.add_argument("--barcode-positions", help="FilterBarcodes CSV. Defaults to <sampledir>/06.segment/mask/<sample>_FilterBarcodes.csv.")
    parser.add_argument("--count-detail", help="05.count/<sample>_count_detail.txt.")
    parser.add_argument("--tissue-bbox", help="CSV containing bbox as x0,y0,width,height.")
    parser.add_argument("--barcode-allowlist", action="append", default=[], help="Optional barcode allowlist; can be used multiple times.")
    parser.add_argument("--barcode-mask", help="Optional tissue mask; auto-detects cavity_filter mask if present.")
    parser.add_argument("--barcode-mask-keep-value", choices=["white", "black"], default="white")
    parser.add_argument("--barcode-mask-threshold", type=int, default=127)
    parser.add_argument("--barcode-in-tissue-only", action="store_true")
    parser.add_argument("--expand-pixels", type=float, default=8.0)
    parser.add_argument("--rounding", choices=["floor", "round"], default="floor")
    parser.add_argument("--chunksize", type=int, default=1_000_000)
    parser.add_argument("--min-umi", type=int, default=10)
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--max-umi", type=int, default=None)
    parser.add_argument("--visualization-preview-dim", type=int, default=1500)
    parser.add_argument("--no-label-output", action="store_true", help="Do not write expanded/QC label TIFF masks.")
