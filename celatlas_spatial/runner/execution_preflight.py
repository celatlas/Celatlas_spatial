"""Preflight checks for real Python step execution."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .executor import StepJobPlan
from .steps import StepPlan


class ExecutionPreflightError(FileNotFoundError):
    """Raised when a selected real-execution step is missing required inputs."""


@dataclass(frozen=True)
class ExecutionPreflightCheck:
    step: str
    label: str
    path: Path
    kind: str


@dataclass(frozen=True)
class ExecutionOutputCollision:
    step: str
    path: Path


@dataclass(frozen=True)
class ExecutionPreflightResult:
    checks: tuple[ExecutionPreflightCheck, ...]
    missing: tuple[ExecutionPreflightCheck, ...]
    output_collisions: tuple[ExecutionOutputCollision, ...] = ()


DIR_OPTIONS = {
    "--genomeDir": "reference genome",
    "--bin-segment-dir": "bin segmentation directory",
    "--runtime-input-dir": "runtime input directory",
    "--square_bin_dir": "square bin directory",
    "--cell-seg-dir": "cell segmentation directory",
    "--matrix_file": "matrix directory",
}

FILE_OPTIONS = {
    "--bam": "featureCounts BAM",
    "--fq": "FASTQ input",
    "--fq1": "FASTQ R1 input",
    "--fq2": "FASTQ R2 input",
    "--count_detail": "count detail",
    "--tif": "image input",
    "--whitelist": "barcode whitelist",
}


def preflight_python_step_jobs(
    step_jobs: list[StepJobPlan],
    *,
    allow_existing_outputs: bool = False,
) -> list[ExecutionPreflightCheck]:
    result = collect_python_step_job_checks(step_jobs)
    errors: list[str] = []
    if result.missing:
        errors.extend(
            f"{check.step}: missing {check.label}: {check.path}"
            for check in result.missing
        )
    if result.output_collisions and not allow_existing_outputs:
        errors.extend(
            f"{collision.step}: existing step output: {collision.path}"
            for collision in result.output_collisions
        )
    if errors:
        raise ExecutionPreflightError("Python step execution preflight failed:\n" + "\n".join(errors))
    return list(result.checks)


def collect_python_step_job_checks(step_jobs: list[StepJobPlan]) -> ExecutionPreflightResult:
    checks: list[ExecutionPreflightCheck] = []
    missing: list[ExecutionPreflightCheck] = []
    output_collisions: list[ExecutionOutputCollision] = []
    for job in step_jobs:
        planned_paths: list[Path] = list(job.planned_existing_paths)
        for step in job.steps:
            decision = job.resume_decisions.get(step.name)
            if decision and decision.action == "skip":
                planned_paths.extend(step.outputs)
                continue
            for output in step.outputs:
                if (
                    output.exists()
                    and not _reanalysis_replaces_outputs(job)
                    and not _is_overwritable_failure_report(step.name, output)
                ):
                    output_collisions.append(ExecutionOutputCollision(step.name, output))
            for check in _checks_for_step(step):
                checks.append(check)
                if _covered_by_planned_path(check.path, planned_paths):
                    continue
                if check.label == "report sample directory" and _contains_planned_path(check.path, planned_paths):
                    continue
                if check.kind == "dir" and not check.path.is_dir():
                    missing.append(check)
                elif check.kind == "file" and not check.path.is_file():
                    missing.append(check)
            planned_paths.extend(step.outputs)
    return ExecutionPreflightResult(tuple(checks), tuple(missing), tuple(output_collisions))


def _reanalysis_replaces_outputs(job: StepJobPlan) -> bool:
    return job.metadata.get("pipeline") == "reanalysis"


def _checks_for_step(step: StepPlan) -> list[ExecutionPreflightCheck]:
    command = step.command
    checks: list[ExecutionPreflightCheck] = []
    for option, label in DIR_OPTIONS.items():
        for value in _option_values(command, option):
            checks.append(ExecutionPreflightCheck(step.name, label, Path(value).expanduser(), "dir"))
    for option, label in FILE_OPTIONS.items():
        for value in _option_values(command, option):
            checks.append(ExecutionPreflightCheck(step.name, label, Path(value).expanduser(), "file"))

    # binSegment takes a local .pth path; stardistCellSegment takes a pretrained model name.
    if step.name == "06.segment/01.binsegment":
        for value in _option_values(command, "--model"):
            checks.append(
                ExecutionPreflightCheck(step.name, "segmentation model", Path(value).expanduser(), "file")
            )

    for value in _option_values(command, "--input"):
        if step.name == "04.featureCounts":
            checks.append(ExecutionPreflightCheck(step.name, "aligned BAM", Path(value).expanduser(), "file"))
        else:
            checks.append(ExecutionPreflightCheck(step.name, "step input directory", Path(value).expanduser(), "dir"))

    if step.name.startswith("06."):
        sample = _first_option_value(command, "--sample")
        input_dir = _first_option_value(command, "--input")
        if sample and input_dir:
            checks.append(
                ExecutionPreflightCheck(
                    step.name,
                    "spatial barcode positions",
                    Path(input_dir).expanduser() / f"{sample}_FilterBarcodes.csv",
                    "file",
                )
            )

    if step.name == "08.report" and command:
        sampledir = _report_sampledir(command)
        if sampledir:
            checks.append(ExecutionPreflightCheck(step.name, "report sample directory", sampledir, "dir"))

    return checks


def _option_values(command: list[str], option: str) -> list[str]:
    values: list[str] = []
    index = 0
    while index < len(command):
        item = command[index]
        if item == option and index + 1 < len(command):
            values.extend(_split_path_list(command[index + 1]))
            index += 2
            continue
        if item.startswith(f"{option}="):
            values.extend(_split_path_list(item.split("=", 1)[1]))
        index += 1
    return values


def _split_path_list(value: str) -> list[str]:
    return [item for item in (part.strip() for part in value.split(",")) if item]


def _first_option_value(command: list[str], option: str) -> str | None:
    values = _option_values(command, option)
    return values[0] if values else None


def _report_sampledir(command: list[str]) -> Path | None:
    executable = Path(command[0]).name
    if executable in {"celatlas_spatial_report", "celatlas_scrna_report"} and len(command) >= 2:
        return Path(command[1]).expanduser()
    return None


def _covered_by_planned_path(path: Path, planned_paths: list[Path]) -> bool:
    resolved = path
    for planned in planned_paths:
        if path == planned:
            return True
        try:
            resolved.relative_to(planned)
            return True
        except ValueError:
            continue
    return False


def _contains_planned_path(path: Path, planned_paths: list[Path]) -> bool:
    for planned in planned_paths:
        try:
            planned.relative_to(path)
            return True
        except ValueError:
            continue
    return False


def _is_overwritable_failure_report(step_name: str, path: Path) -> bool:
    if step_name != "08.report" or not path.is_file():
        return False

    sample = _sample_from_report_name(path.name)
    if not sample:
        return False

    summary_path = path.parent / f"{sample}_pipeline_failure_summary.json"
    if summary_path.is_file() and _summary_points_to_report(summary_path, path):
        return True

    try:
        prefix = path.read_text(encoding="utf-8", errors="ignore")[:8192]
    except OSError:
        return False
    return "Celatlas Pipeline Failure Report" in prefix


def _sample_from_report_name(name: str) -> str | None:
    for suffix in ("_spatial_analysis_report.html", "_scrna_analysis_report.html"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def _summary_points_to_report(summary_path: Path, report_path: Path) -> bool:
    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    output_report = payload.get("output_report")
    if not output_report:
        return False
    return Path(output_report).expanduser() == report_path
