"""Workflow registry for Celatlas runner planning."""

from __future__ import annotations

from dataclasses import dataclass


COMMON_VALUE_FIELDS = {
    "reference_dir",
    "mask_dir",
    "image_dir",
    "thread",
    "featurecounts_thread",
    "featurecounts_threads",
    "bin",
    "cell_num",
    "pixel_size",
    "pixelsize",
    "gem_bin_size",
    "ssdna_threshold_scale",
    "ssdna_mask_expand_pixels",
    "ssdna_min_hole_area",
    "cluster_resolution",
    "cell_cluster_resolution",
    "umi_min_threshold",
    "registration_type",
}
UPSTREAM_VALUE_FIELDS = COMMON_VALUE_FIELDS | {
    "fastq_dir",
    "fastq_name",
    "insert_r2",
    "insertr2",
}
ST_VALUE_FIELDS = {
    "star_match_min",
    "star_match_ratio",
    "star_score_ratio",
    "bbv4_strna_barcode_mismatch",
    "bbv4_strna_mismatch",
}
ST_SX_VALUE_FIELDS = {"star_multimap"}
SX_VALUE_FIELDS = {
    "cell_n_neighbors",
    "cell_n_pcs",
    "cell_min_genes",
    "cell_min_counts",
    "cell_max_genes",
    "cell_max_counts",
    "cell_max_mt",
}
SN_VALUE_FIELDS = {
    "sn_barcode_lownum",
    "sn_cutadapt_min_length",
    "sn_cutadapt_nextseq_trim",
    "sn_cutadapt_overlap",
    "sn_cutadapt_param",
    "sn_star_multimap",
    "sn_strna_star_match_min",
    "sn_strna_star_match_ratio",
    "sn_strna_star_score_ratio",
    "sn_scrna_star_match_min",
    "sn_scrna_star_match_ratio",
    "sn_scrna_star_score_ratio",
    "sn_feature_type",
    "sn_featurecounts_param",
}
REANALYSIS_VALUE_FIELDS = {
    "workspace",
    "genome_dir",
    "genomedir",
    "sampledir",
    "sample_dir",
    "src_dir",
}
POST_STARDIST_VALUE_FIELDS = {
    "cell_segmentation_preset",
    "stardist_python",
    "stardist_labels",
    "stardist_image",
    "stardist_tissue_bbox",
    "stardist_model",
    "stardist_model_dir",
    "stardist_prob_thresh",
    "stardist_nms_thresh",
    "stardist_max_dim",
    "stardist_scale",
    "stardist_n_tiles",
    "stardist_expand_pixels",
    "stardist_min_umi",
    "stardist_min_genes",
    "stardist_tile_size",
    "stardist_tile_overlap",
    "stardist_tile_merge_overlap",
    "stardist_fluorescence_channel",
    "stardist_tile_cache_dir",
    "stardist_exclude_rectangle",
}

ST_SX_FLAG_FIELDS = {"resume_existing"}
ST_FLAG_FIELDS = {"resolve_multigene_umi"}
SN_FLAG_FIELDS = {"sn_resolve_multigene_umi"}
REANALYSIS_FLAG_FIELDS = {"no_clean", "skip_binsegment", "skip_analysis", "skip_report"}
POST_FLAG_FIELDS = {
    "gene_mask_filter",
    "fluorescence_background",
    "enable_cavity_filter",
    "disable_cavity_filter",
    "skip_cavity_filter",
    "cavity_apply",
    "cavity_dry_run",
    "cavity_skip_qc_images",
    "cavity_preserve_gem_support",
    "cavity_force",
    "enable_stardist_cell_segment",
    "enable_stardist",
    "disable_stardist_cell_segment",
    "disable_stardist",
    "stardist_skip_counts",
    "stardist_no_label_output",
    "stardist_tiled_inference",
}


@dataclass(frozen=True)
class WorkflowModule:
    name: str
    backend_script: str
    aliases: tuple[str, ...]
    value_fields: frozenset[str]
    flag_fields: frozenset[str]

    def allows_value_field(self, field: str) -> bool:
        return field in self.value_fields

    def allows_flag_field(self, field: str) -> bool:
        return field in self.flag_fields


WORKFLOW_MODULES = {
    "ST": WorkflowModule(
        name="ST",
        backend_script="Celatlas.sh",
        aliases=("st", "ff", "fresh", "fresh_frozen", "fresh-frozen", "standard", "std", "normal"),
        value_fields=frozenset(UPSTREAM_VALUE_FIELDS | ST_VALUE_FIELDS | ST_SX_VALUE_FIELDS | SX_VALUE_FIELDS | REANALYSIS_VALUE_FIELDS | POST_STARDIST_VALUE_FIELDS),
        flag_fields=frozenset(ST_SX_FLAG_FIELDS | ST_FLAG_FIELDS | POST_FLAG_FIELDS | REANALYSIS_FLAG_FIELDS),
    ),
    "SX": WorkflowModule(
        name="SX",
        backend_script="Celatlas_FFPE.sh",
        aliases=("sx", "ffpe", "probe", "panel", "targeted", "targeted_panel", "targeted-panel"),
        value_fields=frozenset(UPSTREAM_VALUE_FIELDS | ST_SX_VALUE_FIELDS | SX_VALUE_FIELDS | REANALYSIS_VALUE_FIELDS | POST_STARDIST_VALUE_FIELDS),
        flag_fields=frozenset(ST_SX_FLAG_FIELDS | POST_FLAG_FIELDS | REANALYSIS_FLAG_FIELDS),
    ),
    "SN": WorkflowModule(
        name="SN",
        backend_script="Celatlas_SN.sh",
        aliases=("sn", "single_nucleus", "single-nucleus", "snrna", "sn_rna", "sn-rna"),
        value_fields=frozenset(UPSTREAM_VALUE_FIELDS | SN_VALUE_FIELDS | REANALYSIS_VALUE_FIELDS),
        flag_fields=frozenset(SN_FLAG_FIELDS | REANALYSIS_FLAG_FIELDS),
    ),
}

_ALIAS_TO_WORKFLOW = {
    alias: module.name
    for module in WORKFLOW_MODULES.values()
    for alias in module.aliases
}


def get_workflow_module(value: str) -> WorkflowModule:
    key = (value or "").strip()
    workflow_name = _ALIAS_TO_WORKFLOW.get(key.lower(), key)
    try:
        return WORKFLOW_MODULES[workflow_name]
    except KeyError as exc:
        expected = ", ".join(WORKFLOW_MODULES)
        raise KeyError(f"Unknown workflow '{value}'. Expected one of: {expected}.") from exc


def workflow_names() -> tuple[str, ...]:
    return tuple(WORKFLOW_MODULES)
