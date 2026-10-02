#!/usr/bin/env bash

# Failure-report helpers for Celatlas shell entrypoints.
# Source this file; do not execute it directly.

celatlas_failure_report_kind() {
  if [[ "${mode:-}" == "scrna" ]]; then
    printf 'scrna'
  else
    printf 'spatial'
  fi
}

celatlas_failure_report_filename() {
  local kind="$1"
  if [[ "$kind" == "scrna" ]]; then
    printf '%s_scrna_analysis_report.html' "${sample:-sample}"
  else
    printf '%s_spatial_analysis_report.html' "${sample:-sample}"
  fi
}

celatlas_failure_partial_report_filename() {
  local kind="$1"
  if [[ "$kind" == "scrna" ]]; then
    printf '%s_partial_scrna_analysis_report.html' "${sample:-sample}"
  else
    printf '%s_partial_spatial_analysis_report.html' "${sample:-sample}"
  fi
}

celatlas_try_generate_partial_report() {
  local kind="$1"
  local partial_name="$2"

  if [[ "${OPT_SKIP_REPORT:-0}" == "1" ]]; then
    return 1
  fi
  if [[ -z "${sampledir:-}" || -z "${sample:-}" || ! -d "${sampledir:-}" ]]; then
    return 1
  fi

  echo "[Failure Report] Trying best-effort partial ${kind} report: ${partial_name}"

  local helper_dir
  helper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  local repo_dir
  repo_dir="$(cd "${helper_dir}/.." && pwd)"
  local pythonpath="${repo_dir}"
  if [[ -n "${PYTHONPATH:-}" ]]; then
    pythonpath="${pythonpath}:${PYTHONPATH}"
  fi

  if [[ "$kind" == "scrna" ]]; then
    if command -v celatlas_scrna_report >/dev/null 2>&1; then
      celatlas_scrna_report \
        "$sampledir" "$sample" \
        --chemistry "${chemistry:-}" \
        --species "${Species:-}" \
        --output-filename "$partial_name" \
        --fast-mode \
        --verbose
    else
      local generator="${CELATLAS_SCRNA_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/scrna_report_generator.py}"
      [[ -f "$generator" ]] || return 1
      PYTHONPATH="$pythonpath" python3 "$generator" \
        "$sampledir" "$sample" \
        --chemistry "${chemistry:-}" \
        --species "${Species:-}" \
        --output-filename "$partial_name" \
        --fast-mode \
        --verbose
    fi
  else
    local cellsegment_args=()
    if [[ "${OPT_ENABLE_STARDIST_CELL_SEGMENT:-0}" != "1" ]]; then
      cellsegment_args+=(--exclude-cell-segmentation)
    fi
    if command -v celatlas_spatial_report >/dev/null 2>&1; then
      celatlas_spatial_report \
        "$sampledir" "$sample" \
        --chemistry "${chemistry:-}" \
        --species "${Species:-}" \
        --method "${method:-}" \
        --output-filename "$partial_name" \
        --fast-mode \
        --verbose \
        "${cellsegment_args[@]}"
    else
      local generator="${CELATLAS_SPATIAL_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/spatial_report_generator.py}"
      [[ -f "$generator" ]] || return 1
      PYTHONPATH="$pythonpath" python3 "$generator" \
        "$sampledir" "$sample" \
        --chemistry "${chemistry:-}" \
        --species "${Species:-}" \
        --method "${method:-}" \
        --output-filename "$partial_name" \
        --fast-mode \
        --verbose \
        "${cellsegment_args[@]}"
    fi
  fi
}

celatlas_generate_failure_report() {
  local exit_code="$1"
  local line_no="$2"
  local failed_command="$3"
  local display_command="$failed_command"

  if [[ "$display_command" == '"$@"' || "$display_command" == "\$@" ]]; then
    display_command="${CELATLAS_LAST_COMMAND:-$display_command}"
  fi

  if [[ "${CELATLAS_FAILURE_REPORT:-1}" == "0" ]]; then
    return 0
  fi
  if [[ "${CELATLAS_FAILURE_REPORT_DONE:-0}" == "1" ]]; then
    return 0
  fi
  export CELATLAS_FAILURE_REPORT_DONE=1

  trap - ERR
  set +e

  if [[ -z "${sampledir:-}" || -z "${sample:-}" ]]; then
    echo "[Failure Report] sampledir/sample is not initialized; cannot write report."
    return 0
  fi
  mkdir -p "$sampledir"

  local kind
  kind="$(celatlas_failure_report_kind)"
  local output_name
  output_name="$(celatlas_failure_report_filename "$kind")"
  local partial_name
  partial_name="$(celatlas_failure_partial_report_filename "$kind")"
  local partial_path=""

  echo ""
  echo "======================================"
  echo "[Failure Report] Pipeline failed; generating diagnostic report"
  echo "======================================"
  echo "[Failure Report] exit_code=${exit_code} line=${line_no}"
  echo "[Failure Report] failed_command=${display_command}"

  celatlas_try_generate_partial_report "$kind" "$partial_name"
  if [[ -f "${sampledir}/${partial_name}" ]]; then
    partial_path="${sampledir}/${partial_name}"
    echo "[Failure Report] Partial report generated: ${partial_path}"
  else
    echo "[Failure Report] Partial report was not generated; using diagnostic fallback only."
  fi

  local helper_dir
  helper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  local repo_dir
  repo_dir="$(cd "${helper_dir}/.." && pwd)"
  local generator="${repo_dir}/celatlas_spatial/tools/pipeline_failure_report.py"

  if [[ ! -f "$generator" ]]; then
    echo "[Failure Report] fallback generator not found: ${generator}"
    return 0
  fi

  python3 "$generator" \
    --sampledir "$sampledir" \
    --sample "$sample" \
    --mode "${mode:-}" \
    --method "${method:-}" \
    --chemistry "${chemistry:-}" \
    --species "${Species:-}" \
    --workflow "${CELATLAS_WORKFLOW:-${workflow:-}}" \
    --log-file "${log_file:-}" \
    --exit-code "$exit_code" \
    --line-no "$line_no" \
    --failed-command "$display_command" \
    --report-kind "$kind" \
    --output-filename "$output_name" \
    --partial-report "$partial_path"

  echo "[Failure Report] Expected report path: ${sampledir}/${output_name}"
  return 0
}

celatlas_on_pipeline_error() {
  local exit_code="$1"
  local line_no="$2"
  local failed_command="$3"
  celatlas_generate_failure_report "$exit_code" "$line_no" "$failed_command"
  exit "$exit_code"
}

celatlas_on_pipeline_exit() {
  local exit_code="$1"
  local line_no="$2"
  local last_command="$3"
  if [[ "$exit_code" == "0" ]]; then
    return 0
  fi
  celatlas_generate_failure_report "$exit_code" "$line_no" "$last_command"
}
