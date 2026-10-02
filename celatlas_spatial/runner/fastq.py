"""FASTQ discovery helpers for Celatlas runner preflight checks."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import get_config_value, load_profile, merge_configs
from .planner import normalize_job_row


class FastqDetectionError(FileNotFoundError):
    """Raised when FASTQ inputs cannot be detected or validated."""


@dataclass(frozen=True)
class FastqInputs:
    fastq_dir: Path
    fastq_name: str
    layout: str
    r1_files: tuple[Path, ...]
    r2_files: tuple[Path, ...]

    @property
    def sample_fq1(self) -> Path:
        return self.r1_files[0]

    @property
    def fq1_csv(self) -> str:
        return ",".join(str(path) for path in self.r1_files)

    @property
    def fq2_csv(self) -> str:
        return ",".join(str(path) for path in self.r2_files)


def detect_fastq_inputs(
    fastq_dir: str | os.PathLike[str],
    fastq_name: str,
    *,
    validate_files: bool = True,
) -> FastqInputs:
    """Detect FASTQ inputs using the same priority as the shell workflows."""

    root = Path(fastq_dir).expanduser()
    if not root.exists():
        raise FastqDetectionError(f"FASTQ directory not found: {root}")

    multi = _detect_multi_lane(root, fastq_name)
    if multi:
        _validate_pair_counts(multi)
        if validate_files:
            _validate_files(multi)
        return multi

    fold = _detect_fold(root, fastq_name)
    if fold:
        if validate_files:
            _validate_files(fold)
        return fold

    simple = _detect_simple(root, fastq_name)
    if simple:
        if validate_files:
            _validate_files(simple)
        return simple

    raise FastqDetectionError(
        "No FASTQ files found for "
        f"fastq_name={fastq_name!r} in {root}. "
        "Tried multi-lane, fold, and simple naming patterns."
    )


def detect_job_fastqs(
    row: dict[str, str],
    config: dict[str, Any] | None = None,
    *,
    repo_dir: str | os.PathLike[str] | None = None,
    profile_dir: str | os.PathLike[str] | None = None,
    validate_files: bool = True,
) -> FastqInputs | None:
    """Resolve and detect FASTQ inputs for one runner row.

    Returns None for reanalysis/report rows because they start from existing outputs.
    """

    row = normalize_job_row(row)
    if row["pipeline"] != "denovo":
        return None

    merged_config = _merge_row_profiles(row, config or {}, repo_dir, profile_dir)
    fastq_name = row.get("fastq_name") or row.get("chip_number")
    if not fastq_name:
        raise FastqDetectionError("Missing fastq_name and chip_number for FASTQ preflight")

    fastq_dir = row.get("fastq_dir") or get_config_value(merged_config, "fastq_dir")
    if not fastq_dir:
        fastq_root = get_config_value(merged_config, "fastq_root")
        if not fastq_root:
            workspace = get_config_value(merged_config, "workspace") or str(Path.home() / "celatlas_spatial")
            fastq_root = str(Path(workspace) / "fastq")
        chemistry_dir = _chemistry_dir(row.get("chemistry", ""))
        fastq_dir = str(Path(str(fastq_root)) / chemistry_dir)

    return detect_fastq_inputs(fastq_dir, fastq_name, validate_files=validate_files)


def _detect_multi_lane(root: Path, fastq_name: str) -> FastqInputs | None:
    r1_files = sorted(
        list(root.glob(f"{fastq_name}_S*_L*_R1_*.fastq.gz"))
        + list(root.glob(f"{fastq_name}_S*_L*_R1_*.fq.gz"))
    )
    r2_files = sorted(
        list(root.glob(f"{fastq_name}_S*_L*_R2_*.fastq.gz"))
        + list(root.glob(f"{fastq_name}_S*_L*_R2_*.fq.gz"))
    )
    if r1_files and r2_files:
        return FastqInputs(root, fastq_name, "multi_lane", tuple(r1_files), tuple(r2_files))
    return None


def _detect_fold(root: Path, fastq_name: str) -> FastqInputs | None:
    r1_files: list[Path] = []
    r2_files: list[Path] = []
    for fold in ("fold1", "fold2", "fold3", "fold4", "fold5"):
        r1 = root / f"{fastq_name}_{fold}_1.fq.gz"
        r2 = root / f"{fastq_name}_{fold}_2.fq.gz"
        if r1.exists() and r2.exists():
            r1_files.append(r1)
            r2_files.append(r2)
    if r1_files:
        return FastqInputs(root, fastq_name, "fold", tuple(r1_files), tuple(r2_files))
    return None


def _detect_simple(root: Path, fastq_name: str) -> FastqInputs | None:
    r1 = root / f"{fastq_name}_1.fq.gz"
    r2 = root / f"{fastq_name}_2.fq.gz"
    if r1.exists() and r2.exists():
        return FastqInputs(root, fastq_name, "simple", (r1,), (r2,))
    return None


def _validate_pair_counts(inputs: FastqInputs) -> None:
    if len(inputs.r1_files) != len(inputs.r2_files):
        raise FastqDetectionError(
            f"FASTQ R1/R2 file count mismatch for {inputs.fastq_name}: "
            f"{len(inputs.r1_files)} R1 vs {len(inputs.r2_files)} R2"
        )


def _validate_files(inputs: FastqInputs) -> None:
    for path in inputs.r1_files + inputs.r2_files:
        if not path.is_file():
            raise FastqDetectionError(f"FASTQ path is not a file: {path}")
        if path.stat().st_size <= 0:
            raise FastqDetectionError(f"FASTQ file is empty: {path}")
        if not os.access(path, os.R_OK):
            raise FastqDetectionError(f"FASTQ file is not readable: {path}")


def _chemistry_dir(chemistry: str) -> str:
    if chemistry in {"BBV3", "BBV3.1"}:
        return "BBV3.4"
    if chemistry in {"BBV4", "BBV4_L9"}:
        return "BBV4"
    return chemistry


def _merge_row_profiles(
    row: dict[str, str],
    config: dict[str, Any],
    repo_dir: str | os.PathLike[str] | None,
    profile_dir: str | os.PathLike[str] | None,
) -> dict[str, Any]:
    profile_value = row.get("profile", "")
    if not profile_value:
        return config

    repo = Path(repo_dir).expanduser() if repo_dir else Path(__file__).resolve().parents[2]
    base_dir = Path(profile_dir).expanduser() if profile_dir else repo / "configs" / "profiles"
    profiles = []
    for raw_name in profile_value.replace(";", ",").split(","):
        name = raw_name.strip()
        if not name:
            continue
        path = Path(name).expanduser()
        if not path.exists():
            path = base_dir / f"{name}.env" if not path.suffix else base_dir / name
        if path.exists():
            profiles.append(load_profile(path))
    return merge_configs(config, *profiles)
