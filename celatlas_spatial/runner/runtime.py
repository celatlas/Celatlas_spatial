"""Runtime input staging for Celatlas workflow jobs."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .methods import requires_he_image, requires_tissue_tif
from .paths import RuntimePaths, find_he_image, resolve_runtime_paths


class RuntimeStageError(FileNotFoundError):
    """Raised when runtime inputs cannot be staged for a job."""


MANUAL_MASK_TEMPLATES = (
    "manual_mask.png",
    "{sample}_manual_mask.png",
    "{sample}_mask_manual.png",
    "mask_manual.png",
    "he_manual_mask.png",
    "{sample}_he_manual_mask.png",
    "{sample}_he_mask_manual.png",
    "he_mask_manual.png",
)


@dataclass(frozen=True)
class StageAction:
    action: str
    label: str
    source: Path | None = None
    target: Path | None = None

    def describe(self) -> str:
        if self.source and self.target:
            return f"{self.label}: {self.source} -> {self.target}"
        if self.target:
            return f"{self.label}: {self.target}"
        if self.source:
            return f"{self.label}: {self.source}"
        return self.label


@dataclass(frozen=True)
class RuntimeStageResult:
    sample: str
    workflow: str
    mode: str
    target_dir: Path | None
    actions: tuple[StageAction, ...]


def stage_runtime_inputs(
    row: dict[str, str],
    config: dict[str, Any] | None = None,
    *,
    dry_run: bool = False,
) -> RuntimeStageResult:
    """Stage mask/image inputs into the workflow runtime input directory.

    This mirrors the current shell behavior:
    - All spatial workflows use ``<sampledir>/06.segment/mask``.
    - Existing targets are preserved, equivalent to ``cp -n``.
    """

    paths = resolve_runtime_paths(row, config)
    actions: list[StageAction] = []

    if paths.mode == "scrna":
        actions.append(StageAction("skip", "spatial input staging skipped for scRNA mode"))
        return RuntimeStageResult(paths.sample, paths.workflow, paths.mode, None, tuple(actions))

    target_dir = _runtime_input_dir(paths)
    _ensure_dir(target_dir, dry_run=dry_run, actions=actions)

    source_filter_barcodes = paths.mask_dir / f"{paths.sample}_FilterBarcodes.csv"
    target_filter_barcodes = target_dir / source_filter_barcodes.name
    runtime_filter_barcodes_exists = target_filter_barcodes.is_file()

    if paths.pipeline == "reanalysis" and runtime_filter_barcodes_exists:
        actions.append(StageAction("exists", "spatial barcode positions", target=target_filter_barcodes))
    else:
        _copy_required(
            source_filter_barcodes,
            target_dir,
            "spatial barcode positions",
            dry_run=dry_run,
            actions=actions,
        )

    _copy_optional(paths.mask_dir / f"{paths.sample}_tissue_bbox.csv", target_dir, "tissue bbox", dry_run, actions)
    _copy_optional(paths.mask_dir / f"{paths.sample}.barcodeToPos.h5", target_dir, "barcode position h5", dry_run, actions)
    _stage_manual_masks(paths, target_dir, dry_run=dry_run, actions=actions)
    _stage_he_image(paths, target_dir, dry_run=dry_run, actions=actions)
    _stage_tissue_image(paths, target_dir, dry_run=dry_run, actions=actions)

    return RuntimeStageResult(paths.sample, paths.workflow, paths.mode, target_dir, tuple(actions))


def _runtime_input_dir(paths: RuntimePaths) -> Path:
    if paths.runtime_mask_dir is None:
        raise RuntimeStageError(f"Cannot resolve runtime mask directory for sample {paths.sample}")
    return paths.runtime_mask_dir


def _stage_manual_masks(
    paths: RuntimePaths,
    target_dir: Path,
    *,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    for template in MANUAL_MASK_TEMPLATES:
        source = paths.mask_dir / template.format(sample=paths.sample)
        _copy_optional(source, target_dir, "manual mask", dry_run, actions)


def _stage_he_image(
    paths: RuntimePaths,
    target_dir: Path,
    *,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    runtime_he = find_he_image(target_dir, paths.sample)
    if runtime_he:
        actions.append(StageAction("exists", "H&E image", target=runtime_he))
        return

    source_he = find_he_image(paths.image_dir, paths.sample)
    if source_he:
        _copy_existing(source_he, target_dir, "H&E image", dry_run, actions)
        return

    if requires_he_image(paths.method):
        raise RuntimeStageError(
            f"Missing H&E image: {paths.image_dir}/{paths.sample}_he.(tif|tiff|png|jpg|jpeg)"
        )


def _stage_tissue_image(
    paths: RuntimePaths,
    target_dir: Path,
    *,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    source = paths.image_dir / f"{paths.sample}.tif"
    target = target_dir / source.name
    if target.is_file():
        actions.append(StageAction("exists", "tissue image", target=target))
        return

    if requires_tissue_tif(paths.method):
        _copy_required(source, target_dir, "tissue image", dry_run=dry_run, actions=actions)


def _ensure_dir(target_dir: Path, *, dry_run: bool, actions: list[StageAction]) -> None:
    if target_dir.is_dir():
        return
    if target_dir.exists():
        raise RuntimeStageError(f"Runtime input path exists but is not a directory: {target_dir}")
    actions.append(StageAction("mkdir", "runtime input directory", target=target_dir))
    if not dry_run:
        target_dir.mkdir(parents=True, exist_ok=True)


def _copy_required(
    source: Path,
    target_dir: Path,
    label: str,
    *,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    if not source.is_file():
        raise RuntimeStageError(f"Missing {label}: {source}")
    _copy_existing(source, target_dir, label, dry_run, actions)


def _copy_optional(
    source: Path,
    target_dir: Path,
    label: str,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    if source.is_file():
        _copy_existing(source, target_dir, label, dry_run, actions)


def _copy_existing(
    source: Path,
    target_dir: Path,
    label: str,
    dry_run: bool,
    actions: list[StageAction],
) -> None:
    target = target_dir / source.name
    if target.exists():
        actions.append(StageAction("exists", label, source=source, target=target))
        return

    actions.append(StageAction("copy", label, source=source, target=target))
    if not dry_run:
        shutil.copy2(source, target)
