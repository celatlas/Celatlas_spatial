"""Command line interface for the Celatlas Python runner."""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

from .config import ConfigError, get_config_value, load_config
from .execution_preflight import (
    ExecutionPreflightError,
    collect_python_step_job_checks,
    preflight_python_step_jobs,
)
from .executor import StepJobPlan, execute_plans, execute_step_jobs
from .fastq import FastqDetectionError, detect_job_fastqs
from .fastq_manifest import FastqManifestError, FastqManifestResult, prepare_fastq_manifest
from .paths import InputPreflightError, preflight_required_inputs
from .planner import PlanError, build_effective_config, build_job_context, build_job_plan
from .resume import build_resume_plan
from .runtime import RuntimeStageError, RuntimeStageResult, stage_runtime_inputs
from .samples import SampleTableError, is_enabled, read_sample_table
from .step_status import evaluate_step_status
from .steps import build_step_plans, filter_step_plans


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def add_runner_subparsers(subparsers: argparse._SubParsersAction) -> None:
    batch_parser = subparsers.add_parser(
        "batch",
        formatter_class=ArgFormatter,
        help="Run Celatlas workflows from a config file and CSV sample table.",
    )
    _add_batch_args(batch_parser)
    batch_parser.set_defaults(func=_main_batch_from_args)

    run_parser = subparsers.add_parser(
        "run",
        formatter_class=ArgFormatter,
        help="Run one Celatlas workflow through the Python runner.",
    )
    _add_run_args(run_parser)
    run_parser.set_defaults(func=_main_run_from_args)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="celatlas_spatial_runner",
        description="Python runner for Celatlas config/CSV workflows.",
        formatter_class=ArgFormatter,
    )
    subparsers = parser.add_subparsers(dest="command")
    add_runner_subparsers(subparsers)
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 2
    try:
        return int(args.func(args) or 0)
    except (
        ConfigError,
        ExecutionPreflightError,
        SampleTableError,
        PlanError,
        FastqDetectionError,
        FastqManifestError,
        InputPreflightError,
        RuntimeStageError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def batch_main(argv: list[str] | None = None) -> int:
    return main(["batch"] + list(argv if argv is not None else sys.argv[1:]))


def _add_common_runner_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", help="Runner config file: YAML, JSON, or simple env format.")
    parser.add_argument("--repo-dir", help="Celatlas source directory containing Celatlas_run.sh.")
    parser.add_argument("--profile-dir", help="Directory for profile names referenced by CSV/config.")
    parser.add_argument("--log-dir", help="Output directory for runner logs and summary.")
    parser.add_argument(
        "--engine",
        choices=("shell", "python"),
        default="shell",
        help="Execution engine. 'shell' delegates to Celatlas_*.sh; 'python' emits Python step plans for migration.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print commands and write summary without running.")
    parser.add_argument("--preflight", action="store_true", help="Validate FASTQ discovery before execution.")
    parser.add_argument(
        "--fastq-manifest",
        action="store_true",
        help="Check FASTQ input manifest state before execution.",
    )
    parser.add_argument(
        "--reset-fastq-outputs",
        action="store_true",
        help="Reset sequencing-dependent outputs when the FASTQ manifest changed.",
    )
    parser.add_argument(
        "--archive-fastq-outputs",
        action="store_true",
        help="Archive stale sequencing-dependent outputs instead of removing them when --reset-fastq-outputs is used.",
    )
    parser.add_argument(
        "--stage-inputs",
        action="store_true",
        help="Stage spatial mask/image inputs with Python before planning execution.",
    )
    parser.add_argument("--show-steps", action="store_true", help="Print Python-generated internal step commands.")
    parser.add_argument("--show-step-status", action="store_true", help="Print expected step output status.")
    parser.add_argument("--show-resume-plan", action="store_true", help="Print Python resume/skip decisions for steps.")
    parser.add_argument("--from-step", help="Python engine step selector to start from, e.g. 06 or 06.segment/01.binsegment.")
    parser.add_argument("--to-step", help="Python engine step selector to stop at, e.g. 08.report.")
    parser.add_argument("--only-steps", help="Comma-separated Python engine step selector list.")
    parser.add_argument("--show-execution-preflight", action="store_true", help="Print real Python step execution input checks.")
    parser.add_argument(
        "--allow-existing-step-outputs",
        action="store_true",
        help="Allow real Python step execution when selected step output paths already exist.",
    )
    parser.add_argument("--stop-on-error", action="store_true", help="Stop batch execution after the first failure.")


def _add_batch_args(parser: argparse.ArgumentParser) -> None:
    _add_common_runner_args(parser)
    parser.add_argument("--samples", "--jobs", required=True, help="CSV sample table.")
    parser.add_argument(
        "--only",
        help="Comma-separated run_id or chip_number list to run from the CSV.",
    )


def _add_run_args(parser: argparse.ArgumentParser) -> None:
    _add_common_runner_args(parser)
    parser.add_argument("--workflow", required=True, help="ST, SX, or SN. FF/FFPE remain accepted as legacy aliases.")
    parser.add_argument("--pipeline", default="denovo", help="denovo, reanalysis, or report.")
    parser.add_argument("--sample-name")
    parser.add_argument("--tissue")
    parser.add_argument("--chip-number", "--chip", required=True)
    parser.add_argument("--casno", required=True)
    parser.add_argument("--chemistry", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--thread")
    parser.add_argument("--star-thread")
    parser.add_argument("--featurecounts-thread")
    parser.add_argument("--bin")
    parser.add_argument("--pixel-size")
    parser.add_argument("--cluster-resolution")
    parser.add_argument("--cell-cluster-resolution")
    parser.add_argument("--insert-r2")
    parser.add_argument("--cell-num")
    parser.add_argument("--cell-min-genes")
    parser.add_argument("--cell-min-counts")
    parser.add_argument("--cell-max-genes")
    parser.add_argument("--cell-max-counts")
    parser.add_argument("--cell-max-mt")
    parser.add_argument("--fastq-dir")
    parser.add_argument("--fastq-name")
    parser.add_argument("--fastq-root")
    parser.add_argument("--mask-dir")
    parser.add_argument("--image-dir")
    parser.add_argument("--reference-dir")
    parser.add_argument("--genome-dir", "--genomeDir", dest="genome_dir")
    parser.add_argument("--workspace")
    parser.add_argument("--results-root")
    parser.add_argument("--sampledir", "--sample-dir", dest="sampledir")
    parser.add_argument("--src-dir")
    parser.add_argument("--profile")
    parser.add_argument("--extra-args")
    parser.add_argument("--resume-existing", action="store_true")
    parser.add_argument("--gene-mask-filter", action="store_true")
    parser.add_argument("--fluorescence-background", action="store_true")
    parser.add_argument("--enable-cavity-filter", action="store_true")
    parser.add_argument("--cavity-apply", action="store_true")
    parser.add_argument("--cavity-min-component-ratio")
    parser.add_argument("--cavity-min-hole-area")
    parser.add_argument("--cavity-close-radius")
    parser.add_argument("--gem-bin-size")
    parser.add_argument("--ssdna-threshold-scale")
    parser.add_argument("--ssdna-mask-expand-pixels")
    parser.add_argument("--ssdna-min-hole-area")
    parser.add_argument("--cavity-preserve-gem-support", action="store_true")
    parser.add_argument(
        "--cavity-force",
        action="store_true",
        help="Regenerate the generated HE mask and replace an existing cavity-filter backup when applying.",
    )
    parser.add_argument("--enable-stardist-cell-segment", action="store_true")
    parser.add_argument("--cell-segmentation-preset", choices=("dapi",))
    parser.add_argument("--stardist-fluorescence-channel", choices=("gray", "red", "green", "blue"))
    parser.add_argument("--no-clean", action="store_true")
    parser.add_argument("--skip-binsegment", action="store_true")
    parser.add_argument("--skip-analysis", action="store_true")
    parser.add_argument("--skip-report", action="store_true")


def _main_batch_from_args(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    rows = read_sample_table(args.samples)
    selected = _parse_only(args.only)
    _ensure_supported_engine_mode(args)

    shell_plans = []
    step_jobs = []
    for row in rows:
        if not is_enabled(row):
            continue
        if selected and row.get("run_id") not in selected and row.get("chip_number") not in selected:
            continue
        effective_config = build_effective_config(
            row,
            config,
            repo_dir=args.repo_dir,
            profile_dir=args.profile_dir,
        )
        if args.engine == "python":
            job_context = build_job_context(
                row,
                config,
                repo_dir=args.repo_dir,
                profile_dir=args.profile_dir,
            )
            shell_plan = None
            run_id = job_context.run_id
        else:
            shell_plan = build_job_plan(
                row,
                config,
                repo_dir=args.repo_dir,
                profile_dir=args.profile_dir,
            )
            job_context = None
            run_id = shell_plan.run_id
        if args.stage_inputs:
            stage_result = stage_runtime_inputs(row, effective_config, dry_run=args.dry_run)
            _print_stage_result(
                run_id,
                stage_result,
                dry_run=args.dry_run,
            )
        else:
            stage_result = None
        if args.fastq_manifest or args.reset_fastq_outputs:
            _print_fastq_manifest_result(
                run_id,
                prepare_fastq_manifest(
                    row,
                    effective_config,
                    dry_run=args.dry_run,
                    reset_outputs=args.reset_fastq_outputs,
                    archive_outputs=args.archive_fastq_outputs,
                ),
                dry_run=args.dry_run,
                reset_outputs=args.reset_fastq_outputs,
            )
        if args.preflight:
            _print_fastq_preflight(row, effective_config, args)
        if args.show_steps:
            _print_step_plans(row, effective_config, args)
        if args.show_step_status:
            _print_step_status(row, effective_config, args)
        if args.show_resume_plan:
            _print_resume_plan(row, effective_config, args)
        if args.engine == "python":
            assert job_context is not None
            step_jobs.append(_build_step_job(row, effective_config, args, job_context=job_context, stage_result=stage_result))
        else:
            assert shell_plan is not None
            shell_plans.append(shell_plan)

    plans = step_jobs if args.engine == "python" else shell_plans
    if not plans:
        raise PlanError("No enabled sample rows matched the requested selection.")

    if args.engine == "python":
        if args.show_execution_preflight:
            _print_execution_preflight(step_jobs)
        _ensure_python_execution_allowed(step_jobs, args)
        results = execute_step_jobs(
            step_jobs,
            log_dir=args.log_dir or _default_log_dir(),
            dry_run=args.dry_run,
            stop_on_error=args.stop_on_error,
            sample_csv=args.samples,
        )
    else:
        results = execute_plans(
            shell_plans,
            log_dir=args.log_dir or _default_log_dir(),
            dry_run=args.dry_run,
            stop_on_error=args.stop_on_error,
            sample_csv=args.samples,
        )
    return _exit_code_from_results(results)


def _main_run_from_args(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    _ensure_supported_engine_mode(args)
    row = {
        "run_id": args.chip_number,
        "workflow": args.workflow,
        "pipeline": args.pipeline,
        "sample_name": args.sample_name or "",
        "tissue": args.tissue or "",
        "chip_number": args.chip_number,
        "casno": args.casno,
        "chemistry": args.chemistry,
        "species": args.species,
        "method": args.method,
        "mode": args.mode,
        "thread": args.thread or "",
        "star_thread": args.star_thread or "",
        "featurecounts_thread": args.featurecounts_thread or "",
        "bin": args.bin or "",
        "pixel_size": args.pixel_size or "",
        "cluster_resolution": args.cluster_resolution or "",
        "cell_cluster_resolution": args.cell_cluster_resolution or "",
        "insert_r2": args.insert_r2 or "",
        "cell_num": args.cell_num or "",
        "cell_min_genes": args.cell_min_genes or "",
        "cell_min_counts": args.cell_min_counts or "",
        "cell_max_genes": args.cell_max_genes or "",
        "cell_max_counts": args.cell_max_counts or "",
        "cell_max_mt": args.cell_max_mt or "",
        "fastq_dir": args.fastq_dir or "",
        "fastq_name": args.fastq_name or "",
        "fastq_root": args.fastq_root or "",
        "mask_dir": args.mask_dir or "",
        "image_dir": args.image_dir or "",
        "reference_dir": args.reference_dir or "",
        "genome_dir": args.genome_dir or "",
        "workspace": args.workspace or "",
        "results_root": args.results_root or "",
        "sampledir": args.sampledir or "",
        "src_dir": args.src_dir or "",
        "profile": args.profile or "",
        "extra_args": args.extra_args or "",
        "resume_existing": "1" if args.resume_existing else "",
        "gene_mask_filter": "1" if args.gene_mask_filter else "",
        "fluorescence_background": "1" if args.fluorescence_background else "",
        "enable_cavity_filter": "1" if args.enable_cavity_filter else "",
        "cavity_apply": "1" if args.cavity_apply else "",
        "cavity_min_component_ratio": args.cavity_min_component_ratio or "",
        "cavity_min_hole_area": args.cavity_min_hole_area or "",
        "cavity_close_radius": args.cavity_close_radius or "",
        "gem_bin_size": args.gem_bin_size or "",
        "ssdna_threshold_scale": args.ssdna_threshold_scale or "",
        "ssdna_mask_expand_pixels": args.ssdna_mask_expand_pixels or "",
        "ssdna_min_hole_area": args.ssdna_min_hole_area or "",
        "cavity_preserve_gem_support": "1" if args.cavity_preserve_gem_support else "",
        "cavity_force": "1" if args.cavity_force else "",
        "enable_stardist_cell_segment": "1" if args.enable_stardist_cell_segment else "",
        "cell_segmentation_preset": args.cell_segmentation_preset or "",
        "stardist_fluorescence_channel": args.stardist_fluorescence_channel or "",
        "no_clean": "1" if args.no_clean else "",
        "skip_binsegment": "1" if args.skip_binsegment else "",
        "skip_analysis": "1" if args.skip_analysis else "",
        "skip_report": "1" if args.skip_report else "",
    }
    effective_config = build_effective_config(
        row,
        config,
        repo_dir=args.repo_dir,
        profile_dir=args.profile_dir,
    )
    if args.engine == "python":
        job_context = build_job_context(row, config, repo_dir=args.repo_dir, profile_dir=args.profile_dir)
        plan = None
        run_id = job_context.run_id
    else:
        plan = build_job_plan(row, config, repo_dir=args.repo_dir, profile_dir=args.profile_dir)
        job_context = None
        run_id = plan.run_id
    if args.stage_inputs:
        stage_result = stage_runtime_inputs(row, effective_config, dry_run=args.dry_run)
        _print_stage_result(
            run_id,
            stage_result,
            dry_run=args.dry_run,
        )
    else:
        stage_result = None
    if args.fastq_manifest or args.reset_fastq_outputs:
        _print_fastq_manifest_result(
            run_id,
            prepare_fastq_manifest(
                row,
                effective_config,
                dry_run=args.dry_run,
                reset_outputs=args.reset_fastq_outputs,
                archive_outputs=args.archive_fastq_outputs,
            ),
            dry_run=args.dry_run,
            reset_outputs=args.reset_fastq_outputs,
        )
    if args.preflight:
        _print_fastq_preflight(row, effective_config, args)
    if args.show_steps:
        _print_step_plans(row, effective_config, args)
    if args.show_step_status:
        _print_step_status(row, effective_config, args)
    if args.show_resume_plan:
        _print_resume_plan(row, effective_config, args)
    if args.engine == "python":
        assert job_context is not None
        step_jobs = [_build_step_job(row, effective_config, args, job_context=job_context, stage_result=stage_result)]
        if args.show_execution_preflight:
            _print_execution_preflight(step_jobs)
        _ensure_python_execution_allowed(step_jobs, args)
        results = execute_step_jobs(
            step_jobs,
            log_dir=args.log_dir or _default_log_dir(),
            dry_run=args.dry_run,
            stop_on_error=args.stop_on_error,
        )
    else:
        assert plan is not None
        results = execute_plans(
            [plan],
            log_dir=args.log_dir or _default_log_dir(),
            dry_run=args.dry_run,
            stop_on_error=args.stop_on_error,
        )
    return _exit_code_from_results(results)


def _parse_only(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def _exit_code_from_results(results) -> int:
    return 1 if any(result.exit_code != 0 for result in results) else 0


def _ensure_supported_engine_mode(args: argparse.Namespace) -> None:
    if _step_filters_set(args) and args.engine != "python":
        raise PlanError("Step filters require --engine python because shell backends always run their own full workflow.")
    if args.show_execution_preflight and args.engine != "python":
        raise PlanError("--show-execution-preflight requires --engine python.")
    if args.show_execution_preflight and not _step_filters_set(args):
        raise PlanError("--show-execution-preflight requires --from-step/--to-step/--only-steps.")
    if args.engine == "python" and not args.dry_run:
        if not _step_filters_set(args):
            raise PlanError("Real Python step execution currently requires --from-step/--to-step/--only-steps.")


def _build_step_job(
    row: dict[str, str],
    config: dict,
    args: argparse.Namespace,
    *,
    job_context,
    stage_result: RuntimeStageResult | None = None,
) -> StepJobPlan:
    resume_decisions = {}
    if _resume_enabled(row, config):
        resume_plan = build_resume_plan(row, config)
        resume_decisions = {decision.step: decision for decision in resume_plan.decisions}
    return StepJobPlan(
        run_id=job_context.run_id,
        workflow=job_context.workflow,
        chip_number=job_context.chip_number,
        steps=_build_step_plans_for_args(row, config, args, validate_fastq=args.preflight),
        env=job_context.env,
        line_number=job_context.line_number,
        resume_decisions=resume_decisions,
        metadata=job_context.metadata,
        planned_existing_paths=_planned_existing_paths(stage_result) if args.dry_run else (),
    )


def _planned_existing_paths(stage_result: RuntimeStageResult | None) -> tuple[Path, ...]:
    if stage_result is None:
        return ()
    paths: list[Path] = []
    if stage_result.target_dir is not None:
        paths.append(stage_result.target_dir)
    for action in stage_result.actions:
        if action.target is not None and action.action in {"mkdir", "copy", "exists"}:
            paths.append(action.target)
    return tuple(_dedupe_paths(paths))


def _dedupe_paths(paths: Iterable[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        normalized = path.expanduser()
        if normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def _resume_enabled(row: dict[str, str], config: dict) -> bool:
    value = row.get("resume_existing") or get_config_value(config, "resume_existing")
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "enable", "enabled"}


def _print_fastq_preflight(row: dict[str, str], config: dict, args: argparse.Namespace) -> None:
    run_id = row.get("run_id") or row.get("chip_number") or "job"
    paths = preflight_required_inputs(row, config)
    print(
        f"[Celatlas preflight] {run_id}: sampledir={paths.sampledir} "
        f"genome={paths.genome_dir}"
    )
    inputs = detect_job_fastqs(
        row,
        config,
        repo_dir=args.repo_dir,
        profile_dir=args.profile_dir,
        validate_files=True,
    )
    if inputs is None:
        print(f"[Celatlas preflight] {run_id}: FASTQ check skipped for {row.get('pipeline') or 'non-denovo'}")
        return
    print(
        f"[Celatlas preflight] {run_id}: FASTQ layout={inputs.layout} "
        f"pairs={len(inputs.r1_files)} dir={inputs.fastq_dir} name={inputs.fastq_name}"
    )


def _print_step_plans(row: dict[str, str], config: dict, args: argparse.Namespace) -> None:
    run_id = row.get("run_id") or row.get("chip_number") or "job"
    plans = _build_step_plans_for_args(row, config, args, validate_fastq=False)
    if not plans:
        print(f"[Celatlas steps] {run_id}: no generated steps")
        return
    print(f"[Celatlas steps] {run_id}: {len(plans)} generated step command(s)")
    for plan in plans:
        print(f"  [{plan.name}] {plan.command_text()}")


def _print_step_status(row: dict[str, str], config: dict, args: argparse.Namespace) -> None:
    run_id = row.get("run_id") or row.get("chip_number") or "job"
    plans = _build_step_plans_for_args(row, config, args, validate_fastq=False)
    statuses = evaluate_step_status(plans)
    if not statuses:
        print(f"[Celatlas step status] {run_id}: no generated steps")
        return
    print(f"[Celatlas step status] {run_id}: {len(statuses)} generated step status record(s)")
    for status in statuses:
        print(f"  [{status.state}] {status.name}")
        for output in status.outputs:
            marker = "exists" if output.exists else "missing"
            print(f"    - {marker}: {output.path}")


def _print_resume_plan(row: dict[str, str], config: dict, args: argparse.Namespace) -> None:
    plan = build_resume_plan(row, config)
    selected_steps = {step.name for step in _build_step_plans_for_args(row, config, args, validate_fastq=False)}
    print(
        f"[Celatlas resume plan] {plan.run_id}: "
        f"workflow={plan.workflow} resume_enabled={int(plan.resume_enabled)} "
        f"fastq_status={plan.fastq_status}"
    )
    for decision in plan.decisions:
        if selected_steps and decision.step not in selected_steps:
            continue
        print(f"  [{decision.action}] {decision.step}: state={decision.state}; {decision.reason}")


def _build_step_plans_for_args(
    row: dict[str, str],
    config: dict,
    args: argparse.Namespace,
    *,
    validate_fastq: bool,
):
    plans = build_step_plans(row, config, validate_fastq=validate_fastq)
    try:
        return filter_step_plans(
            plans,
            from_step=args.from_step,
            to_step=args.to_step,
            only_steps=args.only_steps,
        )
    except ValueError as exc:
        raise PlanError(str(exc)) from exc


def _step_filters_set(args: argparse.Namespace) -> bool:
    return bool(args.from_step or args.to_step or args.only_steps)


def _print_execution_preflight(step_jobs: list[StepJobPlan]) -> None:
    result = collect_python_step_job_checks(step_jobs)
    print(
        f"[Celatlas execution preflight] checked={len(result.checks)} "
        f"missing={len(result.missing)} existing_outputs={len(result.output_collisions)}"
    )
    missing_keys = {(check.step, check.label, check.path) for check in result.missing}
    for check in result.checks:
        status = "missing" if (check.step, check.label, check.path) in missing_keys else "ok"
        print(f"  [{status}] {check.step}: {check.label}: {check.path}")
    for collision in result.output_collisions:
        print(f"  [existing-output] {collision.step}: {collision.path}")


def _ensure_python_execution_allowed(step_jobs: list[StepJobPlan], args: argparse.Namespace) -> None:
    if args.engine != "python" or args.dry_run:
        return
    empty = [job.run_id for job in step_jobs if not job.steps]
    if empty:
        raise PlanError(f"Real Python step execution selected no steps for job(s): {', '.join(empty)}")
    disallowed = [
        f"{job.run_id}:{step.name}"
        for job in step_jobs
        for step in job.steps
        if not _python_real_step_allowed(step.name)
    ]
    if disallowed:
        raise PlanError(
            "Real Python step execution can only run generated pipeline steps. "
            f"Disallowed selected step(s): {', '.join(disallowed)}"
        )
    checks = preflight_python_step_jobs(
        step_jobs,
        allow_existing_outputs=args.allow_existing_step_outputs,
    )
    print(f"[Celatlas execution preflight] checked {len(checks)} input path(s)")


def _python_real_step_allowed(step_name: str) -> bool:
    return step_name.startswith(
        (
            "00.",
            "01.",
            "02.",
            "03.",
            "04.",
            "05.",
            "06.",
            "07.",
            "08.",
        )
    )


def _print_stage_result(run_id: str, result: RuntimeStageResult, *, dry_run: bool) -> None:
    suffix = " dry-run" if dry_run else ""
    target = result.target_dir if result.target_dir is not None else "N/A"
    print(f"[Celatlas stage{suffix}] {run_id}: target={target}")
    for action in result.actions:
        print(f"  [{action.action}] {action.describe()}")


def _print_fastq_manifest_result(
    run_id: str,
    result: FastqManifestResult,
    *,
    dry_run: bool,
    reset_outputs: bool,
) -> None:
    suffix = " dry-run" if dry_run else ""
    print(
        f"[Celatlas FASTQ manifest{suffix}] {run_id}: "
        f"status={result.status} resume_upstream={int(result.resume_upstream)} "
        f"reset_required={int(result.reset_required)}"
    )
    if result.current_manifest_file is not None:
        print(f"  current={result.current_manifest_file}")
    if result.manifest_file is not None:
        print(f"  previous={result.manifest_file}")
    if not result.reset_actions:
        return
    action_word = "planned" if dry_run or not reset_outputs else "applied"
    print(f"  reset_actions={len(result.reset_actions)} ({action_word})")
    for action in result.reset_actions:
        print(f"  [{action.action}] {action.describe()}")


def _default_log_dir() -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("logs") / "runner" / timestamp


if __name__ == "__main__":
    raise SystemExit(main())
