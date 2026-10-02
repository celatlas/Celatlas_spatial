"""FASTQ input manifest checks for Celatlas workflow jobs."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .config import get_config_value
from .fastq import FastqInputs, detect_job_fastqs
from .paths import resolve_runtime_paths
from .planner import normalize_job_row
from .workflows import get_workflow_module


class FastqManifestError(RuntimeError):
    """Raised when FASTQ manifest state cannot be prepared."""


FASTQ_DEPENDENT_OUTPUTS = (
    "00.sample",
    "01.barcode",
    "02.cutadapt",
    "03.star",
    "04.featureCounts",
    "05.count",
    "06.segment/01.binsegment",
    "06.segment/02.cellsegment",
    "07.outs",
    "06_analysis_wrapper",
    "{sample}_spatial_analysis_report.html",
    "{sample}_scrna_analysis_report.html",
    "{sample}_report.html",
    "{sample}_pipeline_failure_report.html",
    "{sample}_pipeline_failure_summary.json",
)

SN_FASTQ_DEPENDENT_OUTPUTS = FASTQ_DEPENDENT_OUTPUTS


@dataclass(frozen=True)
class FastqManifestEntry:
    read: str
    path: Path
    size: int
    mtime: int

    def line(self) -> str:
        return f"{self.read}\t{self.path}\t{self.size}\t{self.mtime}"


@dataclass(frozen=True)
class ResetAction:
    action: str
    path: Path
    target: Path | None = None

    def describe(self) -> str:
        if self.target is not None:
            return f"{self.path} -> {self.target}"
        return str(self.path)


@dataclass(frozen=True)
class FastqManifestResult:
    sample: str
    workflow: str
    status: str
    manifest_file: Path | None
    current_manifest_file: Path | None
    entries: tuple[FastqManifestEntry, ...]
    resume_upstream: bool
    reset_required: bool
    reset_actions: tuple[ResetAction, ...] = ()


def prepare_fastq_manifest(
    row: dict[str, str],
    config: dict[str, Any] | None = None,
    *,
    dry_run: bool = False,
    reset_outputs: bool = False,
    archive_outputs: bool = False,
) -> FastqManifestResult:
    """Build the current FASTQ manifest and optionally reset stale outputs."""

    config = config or {}
    row = normalize_job_row(row)
    workflow = get_workflow_module(row.get("workflow", ""))
    if row["pipeline"] != "denovo":
        return FastqManifestResult(
            sample=row.get("chip_number") or row.get("sample") or "",
            workflow=workflow.name,
            status="skipped",
            manifest_file=None,
            current_manifest_file=None,
            entries=(),
            resume_upstream=False,
            reset_required=False,
        )

    paths = resolve_runtime_paths(row, config)
    inputs = detect_job_fastqs(row, config, validate_files=True)
    if inputs is None:
        raise FastqManifestError(f"FASTQ inputs unexpectedly skipped for workflow {workflow.name}")

    entries = build_fastq_manifest(inputs)
    manifest_file = paths.sampledir / ".fastq_inputs.tsv"
    current_manifest_file = paths.sampledir / ".fastq_inputs.current.tsv"
    status = fastq_manifest_status(manifest_file, entries)
    has_outputs = has_fastq_dependent_outputs(paths.sampledir, paths.sample, workflow=workflow.name)
    reset_required = (
        status in {"added", "reordered", "changed"}
        or (status == "new" and has_outputs)
        # An explicit reset is also a request for a complete denovo rerun when
        # the FASTQ files are unchanged.  Without this branch a completed run
        # with status=same would retain 00-07 and fail execution preflight.
        or (reset_outputs and has_outputs)
    )
    resume_upstream = status == "same" and _truthy(_get(row, config, "resume_existing", default=False))

    reset_actions: tuple[ResetAction, ...] = ()
    if reset_required:
        reset_reason = "manual" if status == "same" and reset_outputs else status
        reset_actions = tuple(
            plan_fastq_dependent_reset(
                paths.sampledir,
                paths.sample,
                workflow=workflow.name,
                reason=reset_reason,
                archive=archive_outputs,
            )
        )
        if reset_outputs:
            apply_fastq_reset(reset_actions, dry_run=dry_run)

    if not dry_run:
        current_manifest_file.parent.mkdir(parents=True, exist_ok=True)
        write_fastq_manifest(entries, current_manifest_file)

    return FastqManifestResult(
        sample=paths.sample,
        workflow=workflow.name,
        status=status,
        manifest_file=manifest_file,
        current_manifest_file=current_manifest_file,
        entries=entries,
        resume_upstream=resume_upstream,
        reset_required=reset_required,
        reset_actions=reset_actions,
    )


def build_fastq_manifest(inputs: FastqInputs) -> tuple[FastqManifestEntry, ...]:
    entries: list[FastqManifestEntry] = []
    for read, files in (("R1", inputs.r1_files), ("R2", inputs.r2_files)):
        for path in files:
            stat = path.stat()
            entries.append(FastqManifestEntry(read, path, stat.st_size, int(stat.st_mtime)))
    return tuple(entries)


def write_fastq_manifest(entries: tuple[FastqManifestEntry, ...], path: Path) -> None:
    path.write_text("".join(f"{entry.line()}\n" for entry in entries))


def commit_fastq_manifest(sampledir: str | Path) -> Path | None:
    """Promote the current FASTQ manifest after a successful denovo run.

    The previous baseline remains untouched while the pipeline is running.  A
    failed run therefore keeps comparing future attempts with the last
    successful FASTQ set, while the ``.current`` file records what was tried.
    """

    if not sampledir:
        return None
    base = Path(sampledir).expanduser()
    current = base / ".fastq_inputs.current.tsv"
    if not current.is_file():
        return None
    manifest = base / ".fastq_inputs.tsv"
    try:
        current.replace(manifest)
    except OSError as exc:
        raise FastqManifestError(
            f"Could not commit FASTQ manifest {current} -> {manifest}: {exc}"
        ) from exc
    return manifest


def fastq_manifest_status(
    previous_manifest: str | Path,
    current_entries: tuple[FastqManifestEntry, ...],
) -> str:
    previous = Path(previous_manifest)
    current_lines = [entry.line() for entry in current_entries]
    if not previous.is_file():
        return "new"

    previous_lines = previous.read_text().splitlines()
    if current_lines == previous_lines:
        return "same"

    current_set = set(current_lines)
    previous_set = set(previous_lines)
    if previous_set.issubset(current_set):
        return "added" if len(current_set - previous_set) > 0 else "reordered"
    return "changed"


def has_fastq_dependent_outputs(
    sampledir: str | Path,
    sample: str,
    *,
    workflow: str = "ST",
) -> bool:
    return any(path.exists() for path in fastq_dependent_paths(sampledir, sample, workflow=workflow))


def fastq_dependent_paths(
    sampledir: str | Path,
    sample: str,
    *,
    workflow: str = "ST",
) -> tuple[Path, ...]:
    base = Path(sampledir)
    templates = SN_FASTQ_DEPENDENT_OUTPUTS if workflow == "SN" else FASTQ_DEPENDENT_OUTPUTS
    return tuple(base / template.format(sample=sample) for template in templates)


def plan_fastq_dependent_reset(
    sampledir: str | Path,
    sample: str,
    *,
    workflow: str = "ST",
    reason: str = "changed",
    archive: bool = False,
    timestamp: str | None = None,
) -> list[ResetAction]:
    base = Path(sampledir)
    archive_root = base / ".rerun_archive" / f"fastq_{reason}_{timestamp or _timestamp()}"
    actions: list[ResetAction] = []
    for path in fastq_dependent_paths(base, sample, workflow=workflow):
        if not path.exists():
            continue
        if archive:
            relative = path.relative_to(base)
            actions.append(ResetAction("archive", path, _unique_archive_target(archive_root, relative)))
        else:
            actions.append(ResetAction("remove", path))
    return actions


def apply_fastq_reset(actions: tuple[ResetAction, ...] | list[ResetAction], *, dry_run: bool = False) -> None:
    for action in actions:
        if dry_run:
            continue
        if action.action == "archive":
            if action.target is None:
                raise FastqManifestError(f"Archive action missing target for {action.path}")
            action.target.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(action.path), str(action.target))
        elif action.action == "remove":
            if action.path.is_dir():
                shutil.rmtree(action.path)
            elif action.path.exists():
                action.path.unlink()
        else:
            raise FastqManifestError(f"Unknown reset action: {action.action}")


def _unique_archive_target(archive_root: Path, name: str | Path) -> Path:
    relative = Path(name)
    target = archive_root / relative
    index = 1
    while target.exists():
        target = archive_root / relative.parent / f"{relative.name}.{index}"
        index += 1
    return target


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _get(row: dict[str, str], config: dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        if row.get(key) not in (None, ""):
            return row[key]
    value = get_config_value(config, *keys)
    return default if value in (None, "") else value


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "enable", "enabled"}
