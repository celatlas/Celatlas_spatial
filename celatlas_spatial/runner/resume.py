"""Resume planning for Python runner step plans."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .config import get_config_value
from .fastq_manifest import FastqManifestResult, prepare_fastq_manifest
from .planner import normalize_job_row
from .step_status import StepStatus, evaluate_step_status
from .steps import build_step_plans
from .workflows import get_workflow_module


UPSTREAM_STEP_NAMES = {
    "00.sample",
    "01.barcode",
    "02.cutadapt",
    "03.star",
    "04.featureCounts",
    "05.count",
}


@dataclass(frozen=True)
class ResumeDecision:
    step: str
    action: str
    state: str
    reason: str
    status: StepStatus


@dataclass(frozen=True)
class ResumePlan:
    run_id: str
    workflow: str
    resume_enabled: bool
    fastq_status: str
    decisions: tuple[ResumeDecision, ...]
    manifest: FastqManifestResult | None = None


def build_resume_plan(row: dict[str, str], config: dict[str, Any] | None = None) -> ResumePlan:
    config = config or {}
    row = normalize_job_row(row)
    workflow = get_workflow_module(row.get("workflow", ""))
    pipeline = row["pipeline"]
    run_id = row.get("run_id") or row.get("chip_number") or "job"
    resume_enabled = _truthy(_get(row, config, "resume_existing", default=False))
    plans = build_step_plans(row, config, validate_fastq=False)
    statuses = evaluate_step_status(plans)

    manifest: FastqManifestResult | None = None
    fastq_status = "skipped"
    if pipeline == "denovo":
        manifest = prepare_fastq_manifest(row, config, dry_run=True)
        fastq_status = manifest.status

    decisions = tuple(
        _decision_for_step(
            status,
            workflow=workflow.name,
            pipeline=pipeline,
            resume_enabled=resume_enabled,
            fastq_status=fastq_status,
        )
        for status in statuses
    )
    return ResumePlan(run_id, workflow.name, resume_enabled, fastq_status, decisions, manifest)


def _decision_for_step(
    status: StepStatus,
    *,
    workflow: str,
    pipeline: str,
    resume_enabled: bool,
    fastq_status: str,
) -> ResumeDecision:
    if pipeline != "denovo":
        return ResumeDecision(status.name, "run", status.state, f"{pipeline} starts from existing outputs", status)

    if status.name not in UPSTREAM_STEP_NAMES:
        return ResumeDecision(status.name, "run", status.state, "downstream stages are regenerated", status)

    if not resume_enabled:
        return ResumeDecision(status.name, "run", status.state, "resume-existing is disabled", status)

    if fastq_status != "same":
        return ResumeDecision(status.name, "run", status.state, f"FASTQ manifest status is {fastq_status}", status)

    if status.state == "complete":
        return ResumeDecision(status.name, "skip", status.state, "upstream output complete and FASTQ manifest unchanged", status)

    return ResumeDecision(status.name, "run", status.state, f"upstream output state is {status.state}", status)


def _get(row: dict[str, str], config: dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        if row.get(key) not in (None, ""):
            return row[key]
    value = get_config_value(config, *keys)
    return default if value in (None, "") else value


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "enable", "enabled"}
