"""Path resolution and input preflight checks for Celatlas runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import get_config_value
from .methods import backend_method_for, requires_he_image, requires_tissue_tif
from .planner import normalize_job_row, normalize_pipeline, normalize_workflow


class InputPreflightError(FileNotFoundError):
    """Raised when required workflow inputs are missing."""


HE_EXTENSIONS = ("tif", "tiff", "png", "jpg", "jpeg")


@dataclass(frozen=True)
class RuntimePaths:
    workflow: str
    pipeline: str
    sample: str
    casno: str
    chemistry: str
    species: str
    method: str
    mode: str
    workspace: Path
    reference_dir: Path
    genome_dir: Path
    mask_dir: Path
    image_dir: Path
    fastq_dir: Path | None
    fastq_name: str | None
    sampledir: Path
    runtime_mask_dir: Path | None
    src_dir: Path
    rawdata_dir: Path | None


def resolve_runtime_paths(row: dict[str, str], config: dict[str, Any] | None = None) -> RuntimePaths:
    config = config or {}
    row = normalize_job_row(row)
    workflow = normalize_workflow(_required(row, "workflow"))
    pipeline = normalize_pipeline(row.get("pipeline"))
    sample = row.get("chip_number") or row.get("sample") or _required(row, "chip_number")
    casno = _required(row, "casno")
    chemistry = _required(row, "chemistry")
    species = _required(row, "species")
    method = backend_method_for(_required(row, "method"))
    mode = _required(row, "mode")

    workspace = Path(str(_row_or_config(row, config, "workspace") or Path.home() / "celatlas_spatial")).expanduser()
    results_root = Path(str(_row_or_config(row, config, "results_root") or workspace / "results")).expanduser()
    mask_dir = Path(str(_row_or_config(row, config, "mask_dir") or workspace / "ST_mask")).expanduser()
    image_dir = Path(str(_row_or_config(row, config, "image_dir") or workspace / "images")).expanduser()
    src_dir = Path(str(_row_or_config(row, config, "src_dir") or workspace / "src")).expanduser()

    reference_value = _row_or_config(row, config, "reference_dir")
    if workflow == "SX" and not row.get("reference_dir"):
        reference_value = get_config_value(config, "sx_reference_dir", "ffpe_reference_dir", "reference_dir")
    if not reference_value:
        reference_value = workspace / "reference"
    reference_dir = Path(str(reference_value)).expanduser()

    if workflow == "SX":
        genome_dir = reference_dir
    elif pipeline == "reanalysis" and (row.get("genome_dir") or row.get("genomedir")):
        genome_dir = Path(row.get("genome_dir") or row.get("genomedir", "")).expanduser()
    else:
        genome_dir = reference_dir / species

    sampledir_value = row.get("targetdir") or row.get("sampledir") or row.get("sample_dir")
    if sampledir_value:
        sampledir = Path(sampledir_value).expanduser()
    else:
        sampledir = results_root / casno / sample

    runtime_mask_dir = None if mode == "scrna" else sampledir / "06.segment" / "mask"
    rawdata_dir = workspace / "rawdata" / sample if workflow == "SN" else None

    fastq_name = row.get("fastq_name") or sample
    fastq_dir = None
    if pipeline == "denovo":
        fastq_dir_value = row.get("fastq_dir") or get_config_value(config, "fastq_dir")
        if not fastq_dir_value:
            fastq_root = _row_or_config(row, config, "fastq_root") or workspace / "fastq"
            fastq_dir_value = Path(str(fastq_root)) / _chemistry_dir(chemistry)
        fastq_dir = Path(str(fastq_dir_value)).expanduser()

    return RuntimePaths(
        workflow=workflow,
        pipeline=pipeline,
        sample=sample,
        casno=casno,
        chemistry=chemistry,
        species=species,
        method=method,
        mode=mode,
        workspace=workspace,
        reference_dir=reference_dir,
        genome_dir=genome_dir,
        mask_dir=mask_dir,
        image_dir=image_dir,
        fastq_dir=fastq_dir,
        fastq_name=fastq_name,
        sampledir=sampledir,
        runtime_mask_dir=runtime_mask_dir,
        src_dir=src_dir,
        rawdata_dir=rawdata_dir,
    )


def preflight_required_inputs(
    row: dict[str, str],
    config: dict[str, Any] | None = None,
    *,
    check_existing_reanalysis: bool = True,
) -> RuntimePaths:
    paths = resolve_runtime_paths(row, config)
    if paths.pipeline != "report":
        _require_dir(paths.genome_dir, "reference genome")

    if paths.pipeline == "report":
        _require_dir(paths.sampledir, "report sample")
        return paths

    if paths.mode == "scrna":
        return paths

    if paths.pipeline == "reanalysis":
        if check_existing_reanalysis:
            _require_dir(paths.sampledir / "05.count", "reanalysis 05.count")
            _require_file(
                paths.sampledir / "05.count" / f"{paths.sample}_count_detail.txt",
                "reanalysis count_detail",
            )
            if paths.runtime_mask_dir and (paths.runtime_mask_dir / f"{paths.sample}_FilterBarcodes.csv").exists():
                _check_runtime_reanalysis_inputs(paths)
            else:
                _check_source_spatial_inputs(paths)
        return paths

    _check_source_spatial_inputs(paths)
    return paths


def find_he_image(image_dir: Path, sample: str) -> Path | None:
    for ext in HE_EXTENSIONS:
        path = image_dir / f"{sample}_he.{ext}"
        if path.is_file():
            return path
    return None


def _check_source_spatial_inputs(paths: RuntimePaths) -> None:
    _require_file(paths.mask_dir / f"{paths.sample}_FilterBarcodes.csv", "spatial barcode positions")
    if requires_tissue_tif(paths.method):
        _require_file(paths.image_dir / f"{paths.sample}.tif", "tissue image")
    if requires_he_image(paths.method) and not find_he_image(paths.image_dir, paths.sample):
        raise InputPreflightError(
            f"Missing H&E image: {paths.image_dir}/{paths.sample}_he.(tif|tiff|png|jpg|jpeg)"
        )


def _check_runtime_reanalysis_inputs(paths: RuntimePaths) -> None:
    assert paths.runtime_mask_dir is not None
    _require_file(
        paths.runtime_mask_dir / f"{paths.sample}_FilterBarcodes.csv",
        "runtime spatial barcode positions",
    )
    if requires_tissue_tif(paths.method):
        _require_file(paths.runtime_mask_dir / f"{paths.sample}.tif", "runtime tissue image")
    if requires_he_image(paths.method):
        if find_he_image(paths.runtime_mask_dir, paths.sample) or find_he_image(paths.image_dir, paths.sample):
            return
        raise InputPreflightError(
            f"Missing runtime H&E image: {paths.runtime_mask_dir}/{paths.sample}_he.(tif|tiff|png|jpg|jpeg)"
        )


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise InputPreflightError(f"Missing {label}: {path}")


def _require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise InputPreflightError(f"Missing {label} directory: {path}")


def _required(row: dict[str, str], key: str) -> str:
    value = row.get(key, "").strip()
    if not value:
        raise InputPreflightError(f"Missing required value '{key}' for path preflight")
    return value


def _row_or_config(row: dict[str, str], config: dict[str, Any], key: str) -> Any:
    if row.get(key):
        return row[key]
    return get_config_value(config, key)


def _chemistry_dir(chemistry: str) -> str:
    if chemistry in {"BBV3", "BBV3.1"}:
        return "BBV3.4"
    if chemistry in {"BBV4", "BBV4_L9"}:
        return "BBV4"
    return chemistry
