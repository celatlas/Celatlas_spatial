"""Execution and logging for Celatlas runner plans."""

from __future__ import annotations

import codecs
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from .failure_report import FailureReportResult, generate_runner_failure_report
from .fastq_manifest import commit_fastq_manifest
from .planner import JobPlan
from .resume import ResumeDecision
from .steps import StepPlan


@dataclass
class StepJobPlan:
    run_id: str
    workflow: str
    chip_number: str
    steps: list[StepPlan]
    env: dict[str, str]
    line_number: str | None = None
    resume_decisions: dict[str, ResumeDecision] = field(default_factory=dict)
    metadata: dict[str, str] = field(default_factory=dict)
    planned_existing_paths: tuple[Path, ...] = ()

    def command_text(self) -> str:
        lines = []
        for step in self.steps:
            decision = self.resume_decisions.get(step.name)
            action = decision.action.upper() if decision else "RUN"
            reason = f" # {decision.reason}" if decision else ""
            lines.append(f"[{action} {step.name}] {step.command_text()}{reason}")
        return "\n".join(lines)


@dataclass
class JobResult:
    plan: JobPlan | StepJobPlan
    status: str
    exit_code: int
    started_at: str
    ended_at: str
    command_log: Path
    stdout_log: Path | None
    failure_report: FailureReportResult | None = None


def _stream_process_output(process: subprocess.Popen, log_handle) -> None:
    """Forward child output without converting carriage returns to newlines.

    ``tqdm`` refreshes progress bars with ``\\r``.  Reading a text-mode pipe
    iteratively makes Python's universal-newline handling turn those refreshes
    into separate lines.  Consume bytes instead and decode incrementally so
    the terminal can keep the progress bar on one line while the log retains
    the original stream semantics.
    """
    assert process.stdout is not None
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    while True:
        chunk = process.stdout.read(8192)
        if not chunk:
            break
        text = decoder.decode(chunk)
        if text:
            sys.stdout.write(text)
            sys.stdout.flush()
            log_handle.write(text)
            log_handle.flush()

    tail = decoder.decode(b"", final=True)
    if tail:
        sys.stdout.write(tail)
        sys.stdout.flush()
        log_handle.write(tail)
        log_handle.flush()


def execute_plans(
    plans: list[JobPlan],
    *,
    log_dir: str | Path,
    dry_run: bool = False,
    stop_on_error: bool = False,
    sample_csv: str | Path | None = None,
) -> list[JobResult]:
    output_dir = Path(log_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    if sample_csv:
        shutil.copy2(sample_csv, output_dir / "samples.csv.copy")

    results: list[JobResult] = []
    for plan in plans:
        result = _execute_one(plan, output_dir, dry_run=dry_run)
        results.append(result)
        if stop_on_error and result.exit_code != 0:
            break

    _write_summary(results, output_dir / "summary.tsv")
    return results


def execute_step_jobs(
    plans: list[StepJobPlan],
    *,
    log_dir: str | Path,
    dry_run: bool = False,
    stop_on_error: bool = False,
    sample_csv: str | Path | None = None,
) -> list[JobResult]:
    """Write/execute Python step plans."""

    output_dir = Path(log_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    if sample_csv:
        shutil.copy2(sample_csv, output_dir / "samples.csv.copy")

    results: list[JobResult] = []
    for plan in plans:
        result = _dry_run_step_job(plan, output_dir) if dry_run else _execute_step_job(plan, output_dir)
        results.append(result)
        if stop_on_error and result.exit_code != 0:
            break

    _write_summary(results, output_dir / "summary.tsv")
    return results


def _execute_one(plan: JobPlan, output_dir: Path, *, dry_run: bool) -> JobResult:
    started = datetime.now().isoformat(timespec="seconds")
    safe_run_id = _safe_name(plan.run_id)
    command_log = output_dir / f"{safe_run_id}.command.txt"
    stdout_log = None if dry_run else output_dir / f"{safe_run_id}.stdout.log"

    command_log.write_text(plan.command_text() + "\n")

    print(f"[Celatlas runner] {plan.run_id}: {plan.command_text()}")
    if dry_run:
        ended = datetime.now().isoformat(timespec="seconds")
        return JobResult(plan, "DRY_RUN", 0, started, ended, command_log, stdout_log)

    exit_code = 0
    with stdout_log.open("w") as log_handle:
        process = subprocess.Popen(
            plan.command,
            env=plan.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=False,
            bufsize=0,
        )
        _stream_process_output(process, log_handle)
        exit_code = process.wait()

    ended = datetime.now().isoformat(timespec="seconds")
    status = "SUCCESS" if exit_code == 0 else "FAILED"
    failure_report = None
    if exit_code == 0 and plan.metadata.get("pipeline") == "denovo":
        manifest = commit_fastq_manifest(plan.metadata.get("sampledir", ""))
        if manifest is not None:
            print(f"[Celatlas runner] {plan.run_id}: FASTQ manifest committed={manifest}")
    if exit_code != 0:
        failure_report = generate_runner_failure_report(
            metadata=plan.metadata,
            command_text=plan.command_text(),
            exit_code=exit_code,
            stdout_log=stdout_log,
            env=plan.env,
            skip_if_existing=False,
        )
        if failure_report.failure_report:
            print(f"[Celatlas runner] {plan.run_id}: failure report={failure_report.failure_report}")
        elif failure_report.message:
            print(f"[Celatlas runner] {plan.run_id}: failure report {failure_report.status}: {failure_report.message}")
    (output_dir / f"{safe_run_id}.status").write_text(f"{status}\n")
    return JobResult(plan, status, exit_code, started, ended, command_log, stdout_log, failure_report)


def _dry_run_step_job(plan: StepJobPlan, output_dir: Path) -> JobResult:
    started = datetime.now().isoformat(timespec="seconds")
    safe_run_id = _safe_name(plan.run_id)
    command_log = output_dir / f"{safe_run_id}.command.txt"
    command_log.write_text(plan.command_text() + "\n")

    print(f"[Celatlas runner] {plan.run_id}: python step plan ({len(plan.steps)} steps)")
    for step in plan.steps:
        decision = plan.resume_decisions.get(step.name)
        action = decision.action.upper() if decision else "RUN"
        if decision:
            print(f"  [{action}] {step.name}: {decision.reason}")
        else:
            print(f"  [{action}] {step.name}")
        print(f"    {step.command_text()}")

    ended = datetime.now().isoformat(timespec="seconds")
    return JobResult(plan, "DRY_RUN", 0, started, ended, command_log, None)


def _execute_step_job(plan: StepJobPlan, output_dir: Path) -> JobResult:
    started = datetime.now().isoformat(timespec="seconds")
    safe_run_id = _safe_name(plan.run_id)
    command_log = output_dir / f"{safe_run_id}.command.txt"
    stdout_log = output_dir / f"{safe_run_id}.stdout.log"
    command_log.write_text(plan.command_text() + "\n")

    print(f"[Celatlas runner] {plan.run_id}: python step execution ({len(plan.steps)} steps)")
    exit_code = 0
    failed_step = None
    failed_command = ""

    with stdout_log.open("w") as log_handle:
        for step in plan.steps:
            decision = plan.resume_decisions.get(step.name)
            if decision and decision.action == "skip":
                message = f"[Celatlas runner] {plan.run_id}: SKIP {step.name}: {decision.reason}\n"
                sys.stdout.write(message)
                log_handle.write(message)
                continue

            if plan.metadata.get("pipeline") == "reanalysis":
                for output in step.outputs:
                    _remove_reanalysis_output(output, step.name, plan.run_id, log_handle)

            message = f"[Celatlas runner] {plan.run_id}: RUN {step.name}: {step.command_text()}\n"
            sys.stdout.write(message)
            log_handle.write(message)
            process = subprocess.Popen(
                step.command,
                env=plan.env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=False,
                bufsize=0,
            )
            _stream_process_output(process, log_handle)
            exit_code = process.wait()
            if exit_code != 0:
                failed_step = step.name
                failed_command = step.command_text()
                break
            missing_outputs = [path for path in step.outputs if not path.exists()]
            if missing_outputs:
                exit_code = 1
                failed_step = step.name
                failed_command = step.command_text()
                message = (
                    f"[Celatlas runner] {plan.run_id}: output check failed for {step.name}; "
                    f"missing: {', '.join(str(path) for path in missing_outputs)}\n"
                )
                sys.stdout.write(message)
                log_handle.write(message)
                break

    ended = datetime.now().isoformat(timespec="seconds")
    status = "SUCCESS" if exit_code == 0 else "FAILED"
    failure_report = None
    if exit_code == 0 and plan.metadata.get("pipeline") == "denovo":
        manifest = commit_fastq_manifest(plan.metadata.get("sampledir", ""))
        if manifest is not None:
            print(f"[Celatlas runner] {plan.run_id}: FASTQ manifest committed={manifest}")
    if exit_code != 0:
        metadata = dict(plan.metadata)
        if failed_step:
            metadata["failed_step"] = failed_step
        failure_report = generate_runner_failure_report(
            metadata=metadata,
            command_text=failed_command or plan.command_text(),
            exit_code=exit_code,
            stdout_log=stdout_log,
            env=plan.env,
            skip_if_existing=False,
        )
        if failure_report.failure_report:
            print(f"[Celatlas runner] {plan.run_id}: failure report={failure_report.failure_report}")
        elif failure_report.message:
            print(f"[Celatlas runner] {plan.run_id}: failure report {failure_report.status}: {failure_report.message}")

    (output_dir / f"{safe_run_id}.status").write_text(f"{status}\n")
    return JobResult(plan, status, exit_code, started, ended, command_log, stdout_log, failure_report)


def _remove_reanalysis_output(path: Path, step_name: str, run_id: str, log_handle) -> None:
    """Remove only the declared output for a reanalysis step before rerunning it."""
    if not path.exists() and not path.is_symlink():
        return
    inference_dir = path / "stardist_inference"
    tile_cache = inference_dir / "tile_cache"
    if (step_name == "06.segment/02.cellsegment" and path.is_dir() and not path.is_symlink()
            and not inference_dir.is_symlink() and tile_cache.is_dir() and not tile_cache.is_symlink()):
        # Retain only inference tile caches. Matrices, bbox, summaries and QC
        # must be rebuilt, not mistaken for fresh downstream results.
        # Inference validates image/parameter fingerprints before reuse.
        message = f"[Celatlas runner] {run_id}: PRESERVE tile cache: {tile_cache}\n"
        sys.stdout.write(message)
        log_handle.write(message)
        for entry in path.iterdir():
            if entry == inference_dir:
                for artifact in inference_dir.iterdir():
                    if artifact != tile_cache:
                        _remove_reanalysis_output(artifact, "cellsegment artifact", run_id, log_handle)
            else:
                _remove_reanalysis_output(entry, "cellsegment artifact", run_id, log_handle)
        return
    message = f"[Celatlas runner] {run_id}: CLEAN {step_name}: {path}\n"
    sys.stdout.write(message)
    log_handle.write(message)
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def _write_summary(results: list[JobResult], path: Path) -> None:
    with path.open("w") as handle:
        handle.write(
            "run_id\tchip_number\tworkflow\tstatus\tstart_time\tend_time\texit_code\tcommand_log\tstdout_log\tfailure_report\tfailure_summary\n"
        )
        for result in results:
            failure_report = result.failure_report
            handle.write(
                "\t".join(
                    [
                        result.plan.run_id,
                        result.plan.chip_number,
                        result.plan.workflow,
                        result.status,
                        result.started_at,
                        result.ended_at,
                        str(result.exit_code),
                        str(result.command_log),
                        str(result.stdout_log or ""),
                        str(failure_report.failure_report if failure_report and failure_report.failure_report else ""),
                        str(failure_report.summary_json if failure_report and failure_report.summary_json else ""),
                    ]
                )
                + "\n"
            )


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value)
