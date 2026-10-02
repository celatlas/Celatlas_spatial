"""Best-effort failure reports for Python runner executions."""

from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True)
class FailureReportResult:
    status: str
    output_report: Path | None = None
    failure_report: Path | None = None
    summary_json: Path | None = None
    message: str = ""


def generate_runner_failure_report(
    *,
    metadata: Mapping[str, str],
    command_text: str,
    exit_code: int,
    stdout_log: Path | None,
    env: Mapping[str, str] | None = None,
    skip_if_existing: bool = True,
) -> FailureReportResult:
    """Generate a fallback failure report when the shell backend fails."""

    sampledir_value = metadata.get("sampledir") or metadata.get("sample_dir")
    sample = metadata.get("sample") or metadata.get("chip_number")
    if not sampledir_value or not sample:
        return FailureReportResult("skipped", message="sampledir/sample metadata is missing")

    sampledir = Path(sampledir_value).expanduser()
    sampledir.mkdir(parents=True, exist_ok=True)
    existing = sampledir / f"{sample}_pipeline_failure_report.html"
    if skip_if_existing and existing.is_file():
        return FailureReportResult(
            "existing",
            output_report=_standard_output_report(sampledir, sample, metadata.get("mode", "")),
            failure_report=existing,
            summary_json=sampledir / f"{sample}_pipeline_failure_summary.json",
            message="shell failure report already exists",
        )

    generator = _failure_generator()
    if not generator.is_file():
        return FailureReportResult("skipped", message=f"failure report generator not found: {generator}")

    mode = metadata.get("mode", "")
    report_kind = "scrna" if mode == "scrna" else "spatial"
    output_name = _standard_output_report(sampledir, sample, mode).name
    log_file = str(stdout_log) if stdout_log else metadata.get("log_file", "")

    command = [
        sys.executable,
        str(generator),
        "--sampledir",
        str(sampledir),
        "--sample",
        sample,
        "--mode",
        mode,
        "--method",
        metadata.get("method", ""),
        "--chemistry",
        metadata.get("chemistry", ""),
        "--species",
        metadata.get("species", ""),
        "--workflow",
        metadata.get("workflow", ""),
        "--log-file",
        log_file,
        "--exit-code",
        str(exit_code),
        "--line-no",
        "runner",
        "--failed-command",
        command_text,
        "--report-kind",
        report_kind,
        "--output-filename",
        output_name,
    ]

    completed = subprocess.run(
        command,
        env=dict(env or {}),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return FailureReportResult(
            "failed",
            message=f"failure report generator exited {completed.returncode}: {completed.stdout[-1000:]}",
        )

    output_report = sampledir / output_name
    failure_report = sampledir / f"{sample}_pipeline_failure_report.html"
    summary_json = sampledir / f"{sample}_pipeline_failure_summary.json"
    return FailureReportResult(
        "generated",
        output_report=output_report if output_report.exists() else None,
        failure_report=failure_report if failure_report.exists() else None,
        summary_json=summary_json if summary_json.exists() else None,
    )


def _failure_generator() -> Path:
    return Path(__file__).resolve().parents[2] / "celatlas_spatial" / "tools" / "pipeline_failure_report.py"


def _standard_output_report(sampledir: Path, sample: str, mode: str) -> Path:
    if mode == "scrna":
        return sampledir / f"{sample}_scrna_analysis_report.html"
    return sampledir / f"{sample}_spatial_analysis_report.html"
