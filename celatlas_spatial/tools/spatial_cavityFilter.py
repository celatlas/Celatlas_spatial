"""Step 06.segment/01: cavity/tissue mask filtering for square-bin matrices."""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import cv2
from PIL import Image

from celatlas_spatial.tools.cavity_filter import filter_square_bins, generate_binary_mask, write_cavity_overlay_qc


def _split_bins(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _infer_sample(args) -> str:
    if args.sample:
        return args.sample
    if args.sampledir:
        return Path(args.sampledir).resolve().name
    if args.square_bin_dir:
        from celatlas_spatial.tools.cavity_filter.filter_bins import infer_sample_from_square_bin

        return infer_sample_from_square_bin(args.square_bin_dir)
    raise ValueError("Cannot infer sample. Provide --sample or --sampledir.")


def _resolve_paths(args) -> dict:
    sample = _infer_sample(args)
    sampledir = Path(args.sampledir).resolve() if args.sampledir else None
    bin_segment_dir = Path(args.bin_segment_dir).resolve() if args.bin_segment_dir else None
    if bin_segment_dir is None and sampledir is not None:
        bin_segment_dir = sampledir / "06.segment" / "01.binsegment"
    if bin_segment_dir is None and args.square_bin_dir:
        bin_segment_dir = Path(args.square_bin_dir).resolve().parent
    if bin_segment_dir is None:
        raise ValueError(
            "Cannot infer binSegment directory. Provide --sampledir, --bin-segment-dir, or --square-bin-dir."
        )

    square_bin_dir = Path(args.square_bin_dir).resolve() if args.square_bin_dir else bin_segment_dir / "square_bin"
    outdir = Path(args.outdir).resolve() if args.outdir else bin_segment_dir / "cavity_filter"
    mask_path = Path(args.mask).resolve() if args.mask else outdir / "mask" / f"{sample}_tissue_mask.png"
    image_path = Path(args.image).resolve() if args.image else _find_registered_image(bin_segment_dir, sample)

    if args.output_square_bin_dir:
        output_square_bin_dir = Path(args.output_square_bin_dir).resolve()
    elif args.apply:
        output_square_bin_dir = outdir / "square_bin_filtered"
    else:
        output_square_bin_dir = bin_segment_dir / "square_bin_cavity_filtered"

    return {
        "sample": sample,
        "sampledir": sampledir,
        "bin_segment_dir": bin_segment_dir,
        "square_bin_dir": square_bin_dir,
        "outdir": outdir,
        "mask_path": mask_path,
        "image_path": image_path,
        "output_square_bin_dir": output_square_bin_dir,
    }


def _find_registered_image(bin_segment_dir: Path, sample: str) -> Path:
    candidates = [
        bin_segment_dir / "images" / f"{sample}_regist.tif",
        bin_segment_dir / "images" / f"{sample}_regist.tiff",
        bin_segment_dir / "images" / f"{sample}_regist.png",
        bin_segment_dir / "images" / "tissue_hires_image.png",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        "Could not find registered image. Provide --image. Tried: "
        + ", ".join(str(path) for path in candidates)
    )


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _apply_square_bin_replacement(square_bin_dir: Path, filtered_dir: Path, backup_dir: Path, force: bool = False) -> None:
    if backup_dir.exists():
        if not force:
            raise FileExistsError(f"Backup square-bin directory already exists: {backup_dir}")
        shutil.rmtree(backup_dir)
    if square_bin_dir.exists():
        shutil.move(str(square_bin_dir), str(backup_dir))
    shutil.move(str(filtered_dir), str(square_bin_dir))


def _preserve_gem_support(mask_path: Path, source_overlay_dir: Path) -> dict:
    """Restore GEM support only inside enclosed holes of the cavity mask."""
    gem_overlay_path = source_overlay_dir / "3c_gem_heatmap_only.png"
    if not gem_overlay_path.exists():
        raise FileNotFoundError(f"GEM support overlay not found: {gem_overlay_path}")
    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    gem_overlay = cv2.imread(str(gem_overlay_path), cv2.IMREAD_COLOR)
    if mask is None or gem_overlay is None:
        raise ValueError("Could not read cavity mask or GEM support overlay")
    if gem_overlay.shape[:2] != mask.shape[:2]:
        gem_overlay = cv2.resize(
            gem_overlay,
            (mask.shape[1], mask.shape[0]),
            interpolation=cv2.INTER_NEAREST,
        )
    # 3c is white outside the original GEM mask; jet colors, including the
    # lowest-expression blue, represent GEM-supported pixels. Restrict the
    # restoration to enclosed holes so border-connected background remains
    # filtered by the cavity mask.
    gem_gray = cv2.cvtColor(gem_overlay, cv2.COLOR_BGR2GRAY)
    _, gem_support = cv2.threshold(gem_gray, 244, 255, cv2.THRESH_BINARY_INV)
    _, cavity_keep = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
    inverse = cv2.bitwise_not(cavity_keep)
    _, labels, stats, _ = cv2.connectedComponentsWithStats(inverse, connectivity=8)
    border_labels = set()
    height, width = inverse.shape[:2]
    for x in range(width):
        border_labels.add(int(labels[0, x]))
        border_labels.add(int(labels[height - 1, x]))
    for y in range(height):
        border_labels.add(int(labels[y, 0]))
        border_labels.add(int(labels[y, width - 1]))
    enclosed_holes = cv2.threshold(inverse, 254, 255, cv2.THRESH_BINARY)[1]
    for label_id in border_labels:
        if label_id <= 0:
            continue
        enclosed_holes[labels == label_id] = 0
    restore = cv2.bitwise_and(gem_support, enclosed_holes)
    merged = cv2.bitwise_or(cavity_keep, restore)
    Image.fromarray(merged).save(mask_path)
    added = restore
    return {
        "source_overlay": str(gem_overlay_path),
        "cavity_mask_pixels": int(cv2.countNonZero(cavity_keep)),
        "gem_support_pixels": int(cv2.countNonZero(gem_support)),
        "enclosed_hole_pixels": int(cv2.countNonZero(enclosed_holes)),
        "merged_mask_pixels": int(cv2.countNonZero(merged)),
        "added_gem_support_pixels": int(cv2.countNonZero(added)),
    }


def run_cavity_filter(args) -> dict:
    paths = _resolve_paths(args)
    sample = paths["sample"]
    outdir = paths["outdir"]
    mask_path = paths["mask_path"]
    gem_support_summary = None

    print(f"Sample: {sample}")
    print(f"Input square_bin: {paths['square_bin_dir']}")
    print(f"Output square_bin: {paths['output_square_bin_dir']}")
    print(f"Mask: {mask_path}")
    print(f"Registered image: {paths['image_path']}")
    print(f"Dry run: {args.dry_run}")
    print(f"Apply replace: {args.apply}")

    mask_summary = None
    if args.mask and not mask_path.exists():
        raise FileNotFoundError(f"Mask file not found: {mask_path}")
    # A forced rerun must not reuse a mask that may have been modified by an
    # earlier GEM-support pass. Explicit --mask remains authoritative.
    regenerate_mask = not mask_path.exists() or (args.force and not args.mask)
    if regenerate_mask:
        if args.no_generate_mask:
            raise FileNotFoundError(f"Mask file not found and --no-generate-mask was set: {mask_path}")
        print("Generating binary tissue mask..." if not mask_path.exists() else "Regenerating binary tissue mask...")
        mask_summary = generate_binary_mask(
            input_image=paths["image_path"],
            output_mask=mask_path,
            mode=args.mask_mode,
            mask_preset=args.mask_preset,
            max_dim=args.mask_max_dim,
            blur_kernel=args.blur_kernel,
            close_radius=args.close_radius,
            open_radius=args.open_radius,
            min_object_area=args.min_object_area,
            min_hole_area=args.min_hole_area,
            min_component_ratio=args.min_component_ratio,
        )
        _write_json(outdir / "summary" / f"{sample}_mask_summary.json", mask_summary)

    if args.preserve_gem_support:
        source_overlay_dir = paths["bin_segment_dir"] / "images" / "overlays"
        gem_support_summary = _preserve_gem_support(mask_path, source_overlay_dir)
        _write_json(outdir / "summary" / f"{sample}_gem_support_summary.json", gem_support_summary)
        print(
            "Preserved GEM support in cavity mask: "
            f"added {gem_support_summary['added_gem_support_pixels']:,} pixels"
        )

    summaries = filter_square_bins(
        square_bin_dir=paths["square_bin_dir"],
        mask_path=mask_path,
        output_square_bin_dir=paths["output_square_bin_dir"],
        sample=sample,
        bins=_split_bins(args.bins),
        remove_mask_value=args.remove_mask_value,
        mask_threshold=args.mask_threshold,
        dilate=args.dilate,
        barcode_radius=args.barcode_radius,
        copy_extra_files=True,
        dry_run=args.dry_run,
        mask_space=args.mask_space,
        hires_image=args.hires_image,
    )

    for summary in summaries:
        label = summary.get("bin", "unknown")
        if summary.get("status") == "skipped":
            print(f"{label}: skipped ({summary.get('reason')})")
            continue
        print(
            f"{label}: kept {summary['kept_barcodes']:,}/{summary['total_barcodes']:,} "
            f"removed {summary['removed_barcodes']:,} ({summary['removed_fraction'] * 100:.2f}%)"
        )

    overlay_summary = None
    if not args.dry_run and not args.skip_qc_images:
        overlay_dir = outdir / "images" / "overlays"
        source_overlay_dir = paths["bin_segment_dir"] / "images" / "overlays"
        try:
            overlay_summary = write_cavity_overlay_qc(
                source_overlay_dir=source_overlay_dir,
                tissue_mask_path=mask_path,
                output_dir=overlay_dir,
                max_dim=args.qc_max_dim,
            )
            print(f"Cavity QC overlays saved to {overlay_summary['output_dir']}")
        except Exception as exc:
            overlay_summary = {"status": "failed", "error": str(exc)}
            print(f"Cavity QC overlay generation failed: {exc}")

    run_summary = {
        "sample": sample,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_square_bin": str(paths["square_bin_dir"]),
        "output_square_bin": str(paths["output_square_bin_dir"]),
        "mask": str(mask_path),
        "registered_image": str(paths["image_path"]),
        "dry_run": bool(args.dry_run),
        "apply": bool(args.apply),
        "mask_summary": mask_summary,
        "gem_support_summary": gem_support_summary,
        "overlay_images": overlay_summary,
        "bin_summaries": summaries,
    }

    if not args.dry_run:
        summary_dir = outdir / "summary"
        for summary in summaries:
            if summary.get("status") != "completed":
                continue
            _write_json(summary_dir / f"{sample}_{summary['bin']}_filter_summary.json", summary)
        _write_json(summary_dir / f"{sample}_cavity_filter_summary.json", run_summary)

    if args.apply and not args.dry_run:
        backup_dir = paths["bin_segment_dir"] / args.backup_dir_name
        _apply_square_bin_replacement(paths["square_bin_dir"], paths["output_square_bin_dir"], backup_dir, force=args.force)
        run_summary["applied_backup_dir"] = str(backup_dir)
        print(f"Applied filtered square_bin. Backup: {backup_dir}")

    return run_summary


def cavityFilter(args):
    run_cavity_filter(args)


def get_opts_cavityFilter(parser, sub_program):
    parser.add_argument("--sampledir", help="Sample output directory containing 06.segment/01.binsegment.")
    parser.add_argument("--sample", help="Sample name. Defaults to basename of --sampledir.")
    parser.add_argument(
        "--bin-segment-dir",
        help="Path to binSegment output. Defaults to <sampledir>/06.segment/01.binsegment.",
    )
    parser.add_argument("--square-bin-dir", help="Input square_bin directory. Defaults to <bin-segment-dir>/square_bin.")
    parser.add_argument("--image", help="Registered HE image used to generate mask. Auto-detected from <bin-segment-dir>/images if omitted.")
    parser.add_argument("--mask", help="Existing binary mask path. If omitted, a mask is generated.")
    parser.add_argument("--outdir", help="Cavity filter working/output directory. Defaults to <bin-segment-dir>/cavity_filter.")
    parser.add_argument("--output-square-bin-dir", help="Filtered square_bin output. Defaults to <bin-segment-dir>/square_bin_cavity_filtered.")
    parser.add_argument("--dry-run", action="store_true", help="Report filtering counts only; do not write filtered bins.")
    parser.add_argument("--apply", action="store_true", help="Replace input square_bin with filtered output after writing a backup.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Regenerate the generated HE mask and allow replacing an existing backup when --apply is used.",
    )
    parser.add_argument("--backup-dir-name", default="square_bin_before_cavity_filter", help="Backup directory name under <bin-segment-dir> for --apply.")
    parser.add_argument("--bins", default="10,20,50,100,Raw", help="Comma-separated bins to process.")

    parser.add_argument("--no-generate-mask", action="store_true", help="Require --mask instead of generating one.")
    parser.add_argument("--mask-mode", choices=["tissue", "cavity"], default="tissue", help="Generated mask mode. Tissue mode uses white=tissue.")
    parser.add_argument(
        "--mask-preset",
        choices=["default", "adipose"],
        default="default",
        help="Generated mask preset. Use adipose to preserve enclosed light fat/adipocyte regions.",
    )
    parser.add_argument("--mask-max-dim", type=int, default=0, help="Downsample image long edge to this size for mask generation; 0 disables downsampling.")
    parser.add_argument("--blur-kernel", type=int, default=1)
    parser.add_argument("--close-radius", type=int, default=20)
    parser.add_argument("--open-radius", type=int, default=3)
    parser.add_argument("--min-object-area", type=int, default=300)
    parser.add_argument("--min-hole-area", type=int, default=500)
    parser.add_argument("--min-component-ratio", type=float, default=0.01)
    parser.add_argument(
        "--preserve-gem-support",
        action="store_true",
        help="Restore GEM support only inside enclosed holes of the HE/cavity mask. Omit for strict cavity filtering.",
    )

    parser.add_argument("--remove-mask-value", choices=["white", "black"], default="black", help="Mask color to remove. For tissue masks, use black.")
    parser.add_argument("--mask-threshold", type=int, default=127)
    parser.add_argument("--dilate", type=int, default=0)
    parser.add_argument("--barcode-radius", type=int, default=0)
    parser.add_argument("--mask-space", choices=["auto", "hires", "native"], default="auto")
    parser.add_argument("--hires-image", help="Optional hires image used for mask-space alignment.")
    parser.add_argument("--skip-qc-images", action="store_true", help="Do not write cavity mask overlay images.")
    parser.add_argument("--qc-max-dim", type=int, default=4000, help="Maximum width/height for cavity mask overlay images.")
