"""Translate runner config and CSV rows into Celatlas shell commands."""

from __future__ import annotations

import os
import shlex
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .config import config_to_env, get_config_value, load_profile, merge_configs
from .methods import normalize_method
from .workflows import get_workflow_module, workflow_names


class PlanError(ValueError):
    """Raised when a job row cannot be converted to a command."""


PIPELINE_ALIASES = {
    "denovo": "denovo",
    "count": "denovo",
    "full": "denovo",
    "reanalysis": "reanalysis",
    "reanalyze": "reanalysis",
    "reanalyse": "reanalysis",
    "re": "reanalysis",
    "manual": "reanalysis",
    "report": "report",
    "mkreport": "report",
}
LEGACY_REANALYSIS_WORKFLOWS = {"reanalysis", "reanalyze", "reanalyse", "re", "manual"}


VALUE_FIELDS = {
    "reference_dir": "--reference_dir",
    "mask_dir": "--mask_dir",
    "image_dir": "--image_dir",
    "fastq_dir": "--fastq_dir",
    "fastq_name": "--fastq_name",
    "thread": "--thread",
    "bin": "--bin",
    "insert_r2": "--insertR2",
    "insertr2": "--insertR2",
    "cell_num": "--cell_num",
    "pixel_size": "--pixelSize",
    "pixelsize": "--pixelSize",
    "cluster_resolution": "--cluster-resolution",
    "star_multimap": "--star-multimap",
    "star_match_min": "--star-match-min",
    "star_match_ratio": "--star-match-ratio",
    "star_score_ratio": "--star-score-ratio",
    "bbv4_strna_barcode_mismatch": "--bbv4-strna-barcode-mismatch",
    "bbv4_strna_mismatch": "--bbv4-strna-barcode-mismatch",
    "cell_cluster_resolution": "--cell-cluster-resolution",
    "cell_n_neighbors": "--cell-n-neighbors",
    "cell_n_pcs": "--cell-n-pcs",
    "workspace": "--workspace",
    "genome_dir": "--genomeDir",
    "genomedir": "--genomeDir",
    "sampledir": "--sampledir",
    "sample_dir": "--sampledir",
    "src_dir": "--src_dir",
    "cavity_bins": "--cavity-bins",
    "cavity_mask_preset": "--cavity-mask-preset",
    "cavity_min_component_ratio": "--cavity-min-component-ratio",
    "cavity_min_hole_area": "--cavity-min-hole-area",
    "cavity_close_radius": "--cavity-close-radius",
    "gem_bin_size": "--gem-bin-size",
    "ssdna_threshold_scale": "--ssdna-threshold-scale",
    "ssdna_mask_expand_pixels": "--ssdna-mask-expand-pixels",
    "ssdna_min_hole_area": "--ssdna-min-hole-area",
    "stardist_python": "--stardist-python",
    "stardist_labels": "--stardist-labels",
    "stardist_image": "--stardist-image",
    "stardist_tissue_bbox": "--stardist-tissue-bbox",
    "stardist_model": "--stardist-model",
    "stardist_model_dir": "--stardist-model-dir",
    "stardist_prob_thresh": "--stardist-prob-thresh",
    "stardist_nms_thresh": "--stardist-nms-thresh",
    "stardist_max_dim": "--stardist-max-dim",
    "stardist_scale": "--stardist-scale",
    "stardist_n_tiles": "--stardist-n-tiles",
    "stardist_expand_pixels": "--stardist-expand-pixels",
    "stardist_min_umi": "--stardist-min-umi",
    "stardist_min_genes": "--stardist-min-genes",
    "stardist_tile_size": "--stardist-tile-size",
    "stardist_tile_overlap": "--stardist-tile-overlap",
    "stardist_tile_merge_overlap": "--stardist-tile-merge-overlap",
    "stardist_fluorescence_channel": "--stardist-fluorescence-channel",
    "stardist_tile_cache_dir": "--stardist-tile-cache-dir",
    "stardist_exclude_rectangle": "--stardist-exclude-rectangle",
}

FLAG_FIELDS = {
    "resume_existing": "--resume-existing",
    "resolve_multigene_umi": "--resolve-multigene-umi",
    "gene_mask_filter": "--gene-mask-filter",
    "enable_cavity_filter": "--enable-cavity-filter",
    "disable_cavity_filter": "--disable-cavity-filter",
    "skip_cavity_filter": "--skip-cavity-filter",
    "cavity_apply": "--cavity-apply",
    "cavity_dry_run": "--cavity-dry-run",
    "cavity_skip_qc_images": "--cavity-skip-qc-images",
    "cavity_preserve_gem_support": "--cavity-preserve-gem-support",
    "cavity_force": "--cavity-force",
    "enable_stardist_cell_segment": "--enable-stardist-cell-segment",
    "enable_stardist": "--enable-stardist-cell-segment",
    "disable_stardist_cell_segment": "--disable-stardist-cell-segment",
    "disable_stardist": "--disable-stardist-cell-segment",
    "stardist_skip_counts": "--stardist-skip-counts",
    "stardist_no_label_output": "--stardist-no-label-output",
    "stardist_tiled_inference": "--stardist-tiled-inference",
    "no_clean": "--no-clean",
    "skip_binsegment": "--skip-binsegment",
    "skip_analysis": "--skip-analysis",
    "skip_report": "--skip-report",
}

ENV_VALUE_FIELDS = {
    "sn_barcode_lownum": "CELATLAS_SN_BARCODE_LOWNUM",
    "sn_cutadapt_min_length": "CELATLAS_SN_CUTADAPT_MIN_LENGTH",
    "sn_cutadapt_nextseq_trim": "CELATLAS_SN_CUTADAPT_NEXTSEQ_TRIM",
    "sn_cutadapt_overlap": "CELATLAS_SN_CUTADAPT_OVERLAP",
    "sn_cutadapt_param": "CELATLAS_SN_CUTADAPT_PARAM",
    "sn_star_multimap": "CELATLAS_SN_STAR_MULTIMAP",
    "sn_strna_star_match_min": "CELATLAS_SN_STRNA_STAR_MATCH_MIN",
    "sn_strna_star_match_ratio": "CELATLAS_SN_STRNA_STAR_MATCH_RATIO",
    "sn_strna_star_score_ratio": "CELATLAS_SN_STRNA_STAR_SCORE_RATIO",
    "sn_scrna_star_match_min": "CELATLAS_SN_SCRNA_STAR_MATCH_MIN",
    "sn_scrna_star_match_ratio": "CELATLAS_SN_SCRNA_STAR_MATCH_RATIO",
    "sn_scrna_star_score_ratio": "CELATLAS_SN_SCRNA_STAR_SCORE_RATIO",
    "sn_feature_type": "CELATLAS_SN_FEATURE_TYPE",
    "sn_featurecounts_param": "CELATLAS_SN_FEATURECOUNTS_PARAM",
    "sn_resolve_multigene_umi": "CELATLAS_SN_RESOLVE_MULTIGENE_UMI",
    "gem_bin_size": "GEM_BIN_SIZE",
    "umi_min_threshold": "UMI_MIN_THRESHOLD",
    "registration_type": "REGISTRATION_TYPE",
    "fluorescence_background": "CELATLAS_FLUORESCENCE_BACKGROUND",
}

RUNNER_EMPTY_CONFIG = Path(__file__).resolve().parents[2] / "configs" / "runner.empty.env"


@dataclass
class JobContext:
    run_id: str
    workflow: str
    chip_number: str
    env: dict[str, str]
    line_number: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)


@dataclass
class JobPlan:
    run_id: str
    workflow: str
    chip_number: str
    command: list[str]
    env: dict[str, str]
    line_number: str | None = None
    metadata: dict[str, str] = field(default_factory=dict)

    def command_text(self) -> str:
        return shlex.join(self.command)


def normalize_workflow(value: str) -> str:
    try:
        return get_workflow_module(value).name
    except KeyError as exc:
        expected = ", ".join(workflow_names())
        raise PlanError(f"Unknown workflow '{value}'. Expected {expected}.") from exc


def normalize_pipeline(value: str | None) -> str:
    key = (value or "denovo").strip().lower()
    try:
        return PIPELINE_ALIASES[key]
    except KeyError as exc:
        expected = ", ".join(("denovo", "reanalysis", "report"))
        raise PlanError(f"Unknown pipeline '{value}'. Expected {expected}.") from exc


def normalize_job_row(row: dict[str, str]) -> dict[str, str]:
    """Return a row using canonical ST/SX/SN workflow and pipeline fields."""

    normalized = dict(row)
    workflow_value = str(normalized.get("workflow") or "").strip()
    if workflow_value.lower() in LEGACY_REANALYSIS_WORKFLOWS:
        normalized["pipeline"] = normalize_pipeline(normalized.get("pipeline") or "reanalysis")
        workflow_value = _legacy_reanalysis_workflow(normalized)
    else:
        normalized["pipeline"] = normalize_pipeline(normalized.get("pipeline") or "denovo")
    normalized["workflow"] = normalize_workflow(workflow_value)
    return normalized


def _legacy_reanalysis_workflow(row: dict[str, str]) -> str:
    explicit = row.get("source_workflow") or row.get("assay") or row.get("workflow_type")
    if explicit:
        return str(explicit)
    chip = str(row.get("chip_number") or row.get("chipname") or "").strip().upper()
    prefix = chip[:2]
    return prefix if prefix in {"ST", "SX", "SN"} else "ST"


def backend_script_for_workflow(workflow: str) -> str:
    try:
        return get_workflow_module(workflow).backend_script
    except KeyError as exc:
        raise PlanError(f"No backend script configured for workflow '{workflow}'") from exc


def build_job_plan(
    row: dict[str, str],
    base_config: dict[str, Any] | None = None,
    *,
    repo_dir: str | os.PathLike[str] | None = None,
    profile_dir: str | os.PathLike[str] | None = None,
    base_env: dict[str, str] | None = None,
) -> JobPlan:
    row = normalize_job_row(row)
    repo = Path(repo_dir).expanduser() if repo_dir else Path(__file__).resolve().parents[2]
    merged_config = build_effective_config(
        row,
        base_config,
        repo_dir=repo,
        profile_dir=profile_dir,
    )

    workflow = normalize_workflow(_required(row, "workflow"))
    pipeline = normalize_pipeline(row.get("pipeline"))
    if _value(row, merged_config, "cell_segmentation_preset"):
        raise PlanError("cell_segmentation_preset requires the Python engine; use public_cli or --engine python.")
    backend = _resolve_backend(workflow, repo, pipeline=pipeline)
    chip_number = _required(row, "chip_number")
    casno = _required(row, "casno")
    chemistry = _required(row, "chemistry")
    species = _required(row, "species")
    method = _required(row, "method")
    mode = _required(row, "mode")
    sample_name = _value(row, merged_config, "sample_name")

    target_args: list[str] = []
    if sample_name:
        target_args.extend([sample_name, chip_number, casno, chemistry, species, method, mode])
    else:
        target_args.extend([chip_number, casno, chemistry, species, method, mode])

    _append_config_defaults(target_args, row, merged_config)

    for field, option in VALUE_FIELDS.items():
        if not _value_field_allowed(field, workflow):
            continue
        value = _workflow_value(row, merged_config, field, workflow)
        if value not in (None, ""):
            target_args.extend([option, str(value)])

    for field, option in FLAG_FIELDS.items():
        if not _flag_field_allowed(field, workflow):
            continue
        value = _value(row, merged_config, field)
        if _truthy(value):
            target_args.append(option)

    extra_args = _value(row, merged_config, "extra_args")
    if extra_args:
        target_args.extend(shlex.split(str(extra_args)))

    env = _build_job_env(row, merged_config, workflow, chip_number, sample_name, base_env)

    run_id = row.get("run_id") or chip_number
    return JobPlan(
        run_id=run_id,
        workflow=workflow,
        chip_number=chip_number,
        command=["bash", str(backend)] + target_args,
        env=env,
        line_number=row.get("_line_number"),
        metadata=_job_metadata(
            row,
            merged_config,
            workflow=workflow,
            pipeline=pipeline,
            chip_number=chip_number,
            sample_name=sample_name,
            casno=casno,
            chemistry=chemistry,
            species=species,
            method=method,
            mode=mode,
            repo=repo,
        ),
    )


def build_job_context(
    row: dict[str, str],
    base_config: dict[str, Any] | None = None,
    *,
    repo_dir: str | os.PathLike[str] | None = None,
    profile_dir: str | os.PathLike[str] | None = None,
    base_env: dict[str, str] | None = None,
) -> JobContext:
    """Build env and metadata for Python-engine jobs without shell backends."""

    row = normalize_job_row(row)
    repo = Path(repo_dir).expanduser() if repo_dir else Path(__file__).resolve().parents[2]
    merged_config = build_effective_config(
        row,
        base_config,
        repo_dir=repo,
        profile_dir=profile_dir,
    )

    workflow = normalize_workflow(_required(row, "workflow"))
    pipeline = normalize_pipeline(row.get("pipeline"))
    chip_number = _required(row, "chip_number")
    casno = _required(row, "casno")
    chemistry = _required(row, "chemistry")
    species = _required(row, "species")
    method = _required(row, "method")
    mode = _required(row, "mode")
    sample_name = _value(row, merged_config, "sample_name")
    env = _build_job_env(row, merged_config, workflow, chip_number, sample_name, base_env)
    run_id = row.get("run_id") or chip_number
    return JobContext(
        run_id=run_id,
        workflow=workflow,
        chip_number=chip_number,
        env=env,
        line_number=row.get("_line_number"),
        metadata=_job_metadata(
            row,
            merged_config,
            workflow=workflow,
            pipeline=pipeline,
            chip_number=chip_number,
            sample_name=sample_name,
            casno=casno,
            chemistry=chemistry,
            species=species,
            method=method,
            mode=mode,
            repo=repo,
        ),
    )


def build_effective_config(
    row: dict[str, str],
    base_config: dict[str, Any] | None = None,
    *,
    repo_dir: str | os.PathLike[str] | None = None,
    profile_dir: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Merge base config, current workflow section, and job profiles."""

    config = base_config or {}
    repo = Path(repo_dir).expanduser() if repo_dir else Path(__file__).resolve().parents[2]
    profiles = _load_profiles(row.get("profile", ""), repo, profile_dir)
    workflow_config = _workflow_section_config(row, config)
    pipeline_config = _pipeline_section_config(row, config)
    method_config = _method_section_config(row, config)
    return merge_configs(config, workflow_config, pipeline_config, method_config, *profiles)


def _workflow_section_config(row: dict[str, str], config: dict[str, Any]) -> dict[str, Any]:
    workflow_value = row.get("workflow", "")
    if not workflow_value:
        return {}
    try:
        workflow = normalize_job_row(row)["workflow"]
    except PlanError:
        return {}

    section_names: list[str]
    if workflow == "ST":
        section_names = ["ff", "fresh_frozen", "st"]
    elif workflow == "SX":
        section_names = ["ffpe", "panel", "targeted", "sx"]
    elif workflow == "SN":
        section_names = ["snrna", "sn_rna", "sn"]
    else:
        section_names = [workflow.lower()]

    merged: dict[str, Any] = {}
    for section in section_names:
        data = config.get(section)
        if isinstance(data, dict):
            merged = merge_configs(merged, data)
    return merged


def _pipeline_section_config(row: dict[str, str], config: dict[str, Any]) -> dict[str, Any]:
    try:
        pipeline = normalize_job_row(row)["pipeline"]
    except PlanError:
        return {}
    section_names = [pipeline]
    if pipeline == "reanalysis":
        section_names.extend(("reanalyze", "reanalyse"))
    elif pipeline == "report":
        section_names.append("mkreport")
    merged: dict[str, Any] = {}
    for section in section_names:
        data = config.get(section)
        if isinstance(data, dict):
            merged = merge_configs(merged, data)
    return merged


def _method_section_config(row: dict[str, str], config: dict[str, Any]) -> dict[str, Any]:
    method_value = row.get("method", "")
    if not method_value:
        return {}
    try:
        method = normalize_method(method_value)
    except ValueError:
        return {}

    methods = config.get("methods")
    if not isinstance(methods, dict):
        return {}

    section_names = [method, method.lower()]
    if method == "gene_expr":
        section_names.extend(("gene-expression", "gene_expression", "expression"))
    elif method == "ssDNA":
        section_names.extend(("ssdna", "ss_dna", "ss-dna"))
    elif method == "HE":
        section_names.extend(("he", "h&e", "h_and_e", "h-e"))

    merged: dict[str, Any] = {}
    for section in section_names:
        data = methods.get(section)
        if isinstance(data, dict):
            merged = merge_configs(merged, data)
    return merged


def _build_job_env(
    row: dict[str, str],
    config: dict[str, Any],
    workflow: str,
    chip_number: str,
    sample_name: Any,
    base_env: dict[str, str] | None,
) -> dict[str, str]:
    env = dict(base_env or os.environ)
    env.update(config_to_env(config))
    env["CELATLAS_WORKFLOW"] = workflow
    env["CELATLAS_PIPELINE"] = normalize_pipeline(row.get("pipeline"))
    env["CELATLAS_SAMPLE_NAME"] = str(sample_name or "")
    env["CELATLAS_CHIP_NUMBER"] = chip_number
    env["CELATLAS_TISSUE"] = str(row.get("tissue") or "")
    if not _value(row, config, "celatlas_config"):
        env["CELATLAS_CONFIG"] = str(RUNNER_EMPTY_CONFIG)
    for key, value in row.items():
        if key.startswith("celatlas_") and value:
            env[key.upper()] = str(value)
        elif key.upper().startswith("CELATLAS_") and value:
            env[key.upper()] = str(value)
    for field, env_key in ENV_VALUE_FIELDS.items():
        if not (_value_field_allowed(field, workflow) or _flag_field_allowed(field, workflow)):
            continue
        value = _value(row, config, field)
        if value not in (None, ""):
            env[env_key] = str(value)
    _apply_env_runtime_paths(env)
    # subprocess.Popen requires every environment value to be text/bytes.
    # YAML numeric defaults (for example gem_bin_size: 10) otherwise remain
    # integers on this path and fail before the child process starts.
    return {str(key): str(value) for key, value in env.items() if value is not None}


def _job_metadata(
    row: dict[str, str],
    config: dict[str, Any],
    *,
    workflow: str,
    pipeline: str,
    chip_number: str,
    sample_name: Any,
    casno: str,
    chemistry: str,
    species: str,
    method: str,
    mode: str,
    repo: Path,
) -> dict[str, str]:
    metadata = {
        "workflow": workflow,
        "pipeline": pipeline,
        "sample": chip_number,
        "sample_name": str(sample_name or ""),
        "tissue": str(row.get("tissue") or ""),
        "chip_number": chip_number,
        "casno": casno,
        "chemistry": chemistry,
        "species": species,
        "method": method,
        "mode": mode,
        "repo_dir": str(repo),
    }
    try:
        from .paths import resolve_runtime_paths

        paths = resolve_runtime_paths(row, config)
    except Exception:
        return metadata
    metadata.update(
        {
            "sampledir": str(paths.sampledir),
            "genome_dir": str(paths.genome_dir),
            "mask_dir": str(paths.mask_dir),
            "image_dir": str(paths.image_dir),
            "log_file": str(paths.sampledir / "pipeline.log"),
        }
    )
    return metadata


def _resolve_backend(workflow: str, repo_dir: Path, *, pipeline: str = "denovo") -> Path:
    if pipeline == "report":
        raise PlanError("The report pipeline is only available through the Python engine.")
    script_name = "Celatlas_reanalysis.sh" if pipeline == "reanalysis" else backend_script_for_workflow(workflow)
    backend = repo_dir / script_name
    if backend.exists():
        return backend
    path_backend = shutil.which(script_name)
    if path_backend:
        return Path(path_backend)
    raise PlanError(f"Backend script '{script_name}' not found under repo dir or PATH: {repo_dir}")


def _append_config_defaults(
    target_args: list[str],
    row: dict[str, str],
    config: dict[str, Any],
) -> None:
    default_map = {
        "thread": "default_thread",
        "bin": "default_bin",
        "pixel_size": "default_pixel_size",
    }
    for field, config_key in default_map.items():
        if row.get(field):
            continue
        value = get_config_value(config, config_key)
        if value not in (None, ""):
            option = VALUE_FIELDS[field]
            target_args.extend([option, str(value)])


def _value_field_allowed(field: str, workflow: str) -> bool:
    try:
        return get_workflow_module(workflow).allows_value_field(field)
    except KeyError as exc:
        raise PlanError(f"No workflow module configured for '{workflow}'") from exc


def _flag_field_allowed(field: str, workflow: str) -> bool:
    try:
        return get_workflow_module(workflow).allows_flag_field(field)
    except KeyError as exc:
        raise PlanError(f"No workflow module configured for '{workflow}'") from exc


def _load_profiles(
    profile_value: str,
    repo_dir: Path,
    profile_dir: str | os.PathLike[str] | None,
) -> list[dict[str, Any]]:
    if not profile_value:
        return []
    base_dir = Path(profile_dir).expanduser() if profile_dir else repo_dir / "configs" / "profiles"
    profiles: list[dict[str, Any]] = []
    for raw_name in profile_value.replace(";", ",").split(","):
        name = raw_name.strip()
        if not name:
            continue
        path = Path(name).expanduser()
        if not path.exists():
            if not path.suffix:
                path = base_dir / f"{name}.env"
            else:
                path = base_dir / name
        if not path.exists():
            raise PlanError(f"Profile not found: {name} (looked under {base_dir})")
        profiles.append(load_profile(path))
    return profiles


def _required(row: dict[str, str], key: str) -> str:
    value = row.get(key, "").strip()
    if not value:
        line = row.get("_line_number", "?")
        raise PlanError(f"Missing required value '{key}' at CSV line {line}")
    return value


def _value(row: dict[str, str], config: dict[str, Any], key: str) -> Any:
    if key in row and row[key] not in (None, ""):
        return row[key]
    return get_config_value(config, key)


def _workflow_value(
    row: dict[str, str],
    config: dict[str, Any],
    key: str,
    workflow: str,
) -> Any:
    if key == "reference_dir" and workflow == "SX":
        if row.get("reference_dir"):
            return row["reference_dir"]
        return get_config_value(config, "sx_reference_dir", "ffpe_reference_dir", "reference_dir")
    return _value(row, config, key)


def _apply_env_runtime_paths(env: dict[str, str]) -> None:
    """Make env-scoped command-line tools win over a polluted user PATH."""

    env_name = env.get("CELATLAS_ENV_NAME") or "celatlas18"
    env_prefix = _resolve_env_prefix(env, env_name)
    if not env_prefix:
        return

    env["CELATLAS_ENV_PREFIX"] = str(env_prefix)
    if "CELATLAS_CONDA_ROOT" not in env and env_prefix.parent.name == "envs":
        env["CELATLAS_CONDA_ROOT"] = str(env_prefix.parent.parent)

    env_bin = env_prefix / "bin"
    # Some shared tools intentionally live outside the analysis environment.
    # Keep the configured order while preserving the selected env as top priority.
    extra_bin_dirs = [
        Path(value).expanduser()
        for value in env.get("CELATLAS_EXTRA_BIN_DIRS", "").split(os.pathsep)
        if value
    ]
    for extra_bin_dir in reversed(extra_bin_dirs):
        if extra_bin_dir.is_dir():
            _prepend_env_path(env, "PATH", extra_bin_dir)
    if env_bin.is_dir() and _truthy(env.get("CELATLAS_USE_ENV_PATH")):
        _prepend_env_path(env, "PATH", env_bin)

    env_lib = env_prefix / "lib"
    if env_lib.is_dir():
        _prepend_env_path(env, "LD_LIBRARY_PATH", env_lib)

    python_bin = env_bin / "python"
    if not env.get("STARDIST_PYTHON") and python_bin.exists():
        env["STARDIST_PYTHON"] = str(python_bin)

    cache_root = Path(env.get("CELATLAS_CACHE_DIR") or "/tmp/celatlas_spatial_cache").expanduser()
    if not env.get("MPLCONFIGDIR"):
        mpl_dir = cache_root / "matplotlib"
        mpl_dir.mkdir(parents=True, exist_ok=True)
        env["MPLCONFIGDIR"] = str(mpl_dir)
    if not env.get("NUMBA_CACHE_DIR"):
        numba_dir = cache_root / "numba"
        numba_dir.mkdir(parents=True, exist_ok=True)
        env["NUMBA_CACHE_DIR"] = str(numba_dir)


def _resolve_env_prefix(env: dict[str, str], env_name: str) -> Path | None:
    explicit = env.get("CELATLAS_ENV_PREFIX")
    if explicit:
        prefix = Path(explicit).expanduser()
        if prefix.is_dir():
            return prefix

    if env_name == "base":
        base = _existing_dir(env.get("CELATLAS_CONDA_ROOT", ""))
        if base:
            return base

    if env.get("CONDA_DEFAULT_ENV") == env_name:
        active = _existing_dir(env.get("CONDA_PREFIX", ""))
        if active:
            return active

    roots: list[Path] = []
    for key in ("CELATLAS_CONDA_ROOT", "MAMBA_ROOT_PREFIX"):
        root = _existing_dir(env.get(key, ""))
        if root:
            roots.append(root)

    conda_exe = env.get("CONDA_EXE")
    if conda_exe:
        conda_root = Path(conda_exe).expanduser().parent.parent
        if conda_root.is_dir():
            roots.append(conda_root)

    home = env.get("HOME")
    if home:
        home_path = Path(home).expanduser()
        roots.extend(
            [
                home_path / "miniforge3",
                home_path / "mambaforge",
                home_path / "miniconda3",
                home_path / "anaconda3",
            ]
        )

    seen: set[Path] = set()
    for root in roots:
        root = root.expanduser()
        if root in seen:
            continue
        seen.add(root)
        prefix = root / "envs" / env_name
        if prefix.is_dir():
            return prefix

    return None


def _existing_dir(value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value).expanduser()
    if path.is_dir():
        return path
    return None


def _prepend_env_path(env: dict[str, str], key: str, path: Path) -> None:
    entry = str(path)
    current = env.get(key, "")
    parts = [part for part in current.split(os.pathsep) if part and part != entry]
    env[key] = os.pathsep.join([entry] + parts)


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on", "enable", "enabled"}
