"""Space Ranger-style public CLI for single-sample Celatlas runs."""

from __future__ import annotations

import argparse
from pathlib import Path

from .cli import ArgFormatter, main as runner_main


PIPELINE_BY_COMMAND = {
    "count": "denovo",
    "reanalyze": "reanalysis",
    "mkreport": "report",
}


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 2
    _validate_args(parser, args)
    return runner_main(_runner_argv(args))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="celatlas",
        description="Celatlas single-sample spatial transcriptomics pipelines.",
        formatter_class=ArgFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")
    for command, help_text in (
        ("count", "Run a complete FASTQ-to-report analysis."),
        ("reanalyze", "Re-run downstream analysis from existing count outputs."),
        ("mkreport", "Generate only the final report from existing outputs."),
    ):
        subparser = subparsers.add_parser(command, help=help_text, formatter_class=ArgFormatter)
        _add_single_sample_args(subparser, command)
    return parser


def _add_single_sample_args(parser: argparse.ArgumentParser, command: str) -> None:
    parser.add_argument("--id", required=True, help="Unique chip/run identifier.")
    parser.add_argument("--sample-name", required=True, help="Biological sample name used in metadata.")
    parser.add_argument("--tissue", required=True, help="Tissue name stored in run and report metadata.")
    parser.add_argument("--targetdir", required=True, help="Absolute final sample output directory.")
    parser.add_argument("--chemistry", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--workflow", required=True, choices=("ST", "SX", "SN"))
    parser.add_argument("--image", required=True, choices=("gene_expr", "ssDNA", "HE"))
    parser.add_argument("--mode", choices=("strna", "scrna"), default="strna")
    parser.add_argument("--project", help="Project/case metadata; defaults to the target directory parent name.")
    parser.add_argument("--config", help="Runner YAML, JSON, or env config.")
    parser.add_argument("--profile")
    parser.add_argument("--fastqs", help="FASTQ directory. Required by count.")
    parser.add_argument("--fastq-name", help="FASTQ filename prefix; defaults to --id.")
    parser.add_argument("--mask-dir")
    parser.add_argument("--image-dir")
    parser.add_argument("--reference-dir")
    parser.add_argument("--thread")
    parser.add_argument("--star-thread", help="STAR alignment thread count; defaults to runner config or --thread.")
    parser.add_argument("--featurecounts-thread")
    parser.add_argument("--bin")
    parser.add_argument("--pixel-size")
    parser.add_argument("--cluster-resolution", help="Leiden resolution for binned spatial clustering.")
    parser.add_argument("--cell-cluster-resolution", help="Leiden resolution for cell clustering.")
    parser.add_argument("--gem-bin-size", help="GEM heatmap aggregation size in microns; 10 gives finer cells than 20.")
    parser.add_argument(
        "--ssdna-threshold-scale",
        help="Scale ssDNA Otsu threshold; below 1.0 retains dim tissue edges.",
    )
    parser.add_argument(
        "--ssdna-mask-expand-pixels",
        help="Expand the final ssDNA tissue mask by this many registered-image pixels.",
    )
    parser.add_argument("--ssdna-min-hole-area", help="Fill only smaller enclosed ssDNA mask holes (pixels).")
    parser.add_argument("--cell-num")
    parser.add_argument("--cell-min-genes")
    parser.add_argument("--cell-min-counts")
    parser.add_argument("--cell-max-genes")
    parser.add_argument("--cell-max-counts")
    parser.add_argument("--cell-max-mt")
    parser.add_argument("--gene-mask-filter", action="store_true")
    parser.add_argument("--fluorescence-background", action="store_true",
                        help="Keep fluorescence tissue backgrounds black when using HE ROI registration.")
    parser.add_argument("--enable-cell-segmentation", action="store_true")
    parser.add_argument("--cell-segmentation-preset", choices=("dapi",),
                        help="HE ROI registration + black background + full-resolution tiled DAPI segmentation and cell analysis.")
    parser.add_argument("--stardist-fluorescence-channel", choices=("gray", "red", "green", "blue"),
                        help="Channel for tiled fluorescence segmentation; DAPI preset defaults to blue (grayscale images also supported).")
    parser.add_argument("--cavity-mode", choices=("off", "qc", "apply"), default="off")
    parser.add_argument("--cavity-min-component-ratio", help="Minimum tissue-component area relative to the largest component.")
    parser.add_argument("--cavity-min-hole-area", help="Fill tissue-mask holes smaller than this area in pixels.")
    parser.add_argument("--cavity-close-radius", help="Tissue-mask closing radius in full-resolution pixels.")
    parser.add_argument(
        "--preserve-gem-support",
        action="store_true",
        help="Restore GEM support only inside enclosed HE-mask holes; omit for strict cavity filtering.",
    )
    parser.add_argument(
        "--cavity-force",
        action="store_true",
        help="Replace the existing cavity-filter backup when applying.",
    )
    parser.add_argument("--log-dir")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--resume-existing",
        action="store_true",
        help="Reuse complete upstream outputs when the FASTQ manifest is unchanged.",
    )
    parser.add_argument("--skip-preflight", action="store_true")
    parser.add_argument("--allow-existing-step-outputs", action="store_true")
    if command == "count":
        parser.add_argument(
            "--reset-fastq-outputs",
            action="store_true",
            help="Run denovo from FASTQ again, resetting sequencing-dependent outputs.",
        )
        parser.add_argument(
            "--archive-fastq-outputs",
            action="store_true",
            help="Archive existing outputs under .rerun_archive before the denovo rerun.",
        )
    if command == "reanalyze":
        parser.add_argument("--skip-binsegment", action="store_true")
        parser.add_argument("--skip-analysis", action="store_true")


def _validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    targetdir = Path(args.targetdir).expanduser()
    if not targetdir.is_absolute():
        parser.error("--targetdir must be an absolute path")
    if args.command == "count" and not args.fastqs:
        parser.error("celatlas count requires --fastqs")
    if args.command == "count" and args.archive_fastq_outputs and not args.reset_fastq_outputs:
        parser.error("--archive-fastq-outputs requires --reset-fastq-outputs")


def _runner_argv(args: argparse.Namespace) -> list[str]:
    targetdir = Path(args.targetdir).expanduser()
    pipeline = PIPELINE_BY_COMMAND[args.command]
    command = [
        "run",
        "--engine",
        "python",
        "--pipeline",
        pipeline,
        "--workflow",
        args.workflow,
        "--chip-number",
        args.id,
        "--sample-name",
        args.sample_name,
        "--tissue",
        args.tissue,
        "--casno",
        args.project or targetdir.parent.name or "celatlas",
        "--chemistry",
        args.chemistry,
        "--species",
        args.species,
        "--method",
        args.image,
        "--mode",
        args.mode,
        "--sampledir",
        str(targetdir),
    ]
    _append_value(command, "--config", args.config)
    _append_value(command, "--profile", args.profile)
    _append_value(command, "--fastq-dir", args.fastqs)
    _append_value(command, "--fastq-name", args.fastq_name or args.id)
    _append_value(command, "--mask-dir", args.mask_dir)
    _append_value(command, "--image-dir", args.image_dir)
    _append_value(command, "--reference-dir", args.reference_dir)
    _append_value(command, "--thread", args.thread)
    _append_value(command, "--star-thread", args.star_thread)
    _append_value(command, "--featurecounts-thread", args.featurecounts_thread)
    _append_value(command, "--bin", args.bin)
    _append_value(command, "--pixel-size", args.pixel_size)
    _append_value(command, "--cluster-resolution", args.cluster_resolution)
    _append_value(command, "--cell-cluster-resolution", args.cell_cluster_resolution)
    _append_value(command, "--gem-bin-size", args.gem_bin_size)
    _append_value(command, "--ssdna-threshold-scale", args.ssdna_threshold_scale)
    _append_value(command, "--ssdna-mask-expand-pixels", args.ssdna_mask_expand_pixels)
    _append_value(command, "--ssdna-min-hole-area", args.ssdna_min_hole_area)
    _append_value(command, "--cell-num", args.cell_num)
    _append_value(command, "--cell-min-genes", args.cell_min_genes)
    _append_value(command, "--cell-min-counts", args.cell_min_counts)
    _append_value(command, "--cell-max-genes", args.cell_max_genes)
    _append_value(command, "--cell-max-counts", args.cell_max_counts)
    _append_value(command, "--cell-max-mt", args.cell_max_mt)
    _append_value(command, "--cavity-min-component-ratio", args.cavity_min_component_ratio)
    _append_value(command, "--cavity-min-hole-area", args.cavity_min_hole_area)
    _append_value(command, "--cavity-close-radius", args.cavity_close_radius)
    _append_value(command, "--log-dir", args.log_dir)
    _append_value(command, "--cell-segmentation-preset", args.cell_segmentation_preset)
    _append_value(command, "--stardist-fluorescence-channel", args.stardist_fluorescence_channel)

    if pipeline == "denovo":
        command.extend(("--from-step", "00", "--to-step", "08.report", "--fastq-manifest", "--stage-inputs"))
    elif pipeline == "reanalysis":
        starts_with_post_binsegment = (
            args.cavity_mode in {"qc", "apply"}
            or args.enable_cell_segmentation
            or args.cell_segmentation_preset
        )
        from_step = "07" if args.skip_binsegment and not starts_with_post_binsegment else "06"
        command.extend(("--from-step", from_step, "--to-step", "08.report", "--stage-inputs"))
    else:
        command.extend(("--only-steps", "08.report", "--allow-existing-step-outputs"))

    if not args.skip_preflight:
        command.append("--preflight")
    if args.dry_run:
        command.extend(("--dry-run", "--show-steps"))
    if args.resume_existing:
        command.append("--resume-existing")
    if getattr(args, "reset_fastq_outputs", False):
        command.append("--reset-fastq-outputs")
    if getattr(args, "archive_fastq_outputs", False):
        command.append("--archive-fastq-outputs")
    if args.allow_existing_step_outputs and "--allow-existing-step-outputs" not in command:
        command.append("--allow-existing-step-outputs")
    if args.gene_mask_filter:
        command.append("--gene-mask-filter")
    if args.fluorescence_background:
        command.append("--fluorescence-background")
    if args.enable_cell_segmentation:
        command.append("--enable-stardist-cell-segment")
    if args.cavity_mode in {"qc", "apply"}:
        command.append("--enable-cavity-filter")
    if args.cavity_mode == "apply":
        command.append("--cavity-apply")
    if args.preserve_gem_support:
        command.append("--cavity-preserve-gem-support")
    if args.cavity_force:
        command.append("--cavity-force")
    if getattr(args, "skip_binsegment", False):
        command.append("--skip-binsegment")
    if getattr(args, "skip_analysis", False):
        command.append("--skip-analysis")
    return command


def _append_value(command: list[str], option: str, value: str | None) -> None:
    if value not in (None, ""):
        command.extend((option, str(value)))


if __name__ == "__main__":
    raise SystemExit(main())
