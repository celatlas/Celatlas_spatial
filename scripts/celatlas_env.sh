#!/usr/bin/env bash

# Shared environment interface for Celatlas v1.8 shell entrypoints.
# Source this file; do not execute it directly.

CELATLAS_REPO_DIR="${CELATLAS_REPO_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
source "${CELATLAS_REPO_DIR}/scripts/celatlas_config.sh"
celatlas_config_load

celatlas_env_try_conda_root() {
  local candidate="${1:-}"
  if [[ -z "${CELATLAS_CONDA_ROOT:-}" && -n "$candidate" && -d "$candidate" ]]; then
    CELATLAS_CONDA_ROOT="$candidate"
  fi
}

celatlas_env_init_defaults() {
  CELATLAS_ENV_NAME="${CELATLAS_ENV_NAME:-celatlas18}"
  CELATLAS_CONDA_ROOT="${CELATLAS_CONDA_ROOT:-}"

  if [[ -z "$CELATLAS_CONDA_ROOT" ]]; then
    if [[ -n "${CONDA_EXE:-}" ]]; then
      celatlas_env_try_conda_root "$(cd "$(dirname "${CONDA_EXE}")/.." && pwd)"
    fi
    celatlas_env_try_conda_root "${MAMBA_ROOT_PREFIX:-}"
    if [[ -z "$CELATLAS_CONDA_ROOT" ]] && command -v conda >/dev/null 2>&1; then
      celatlas_env_try_conda_root "$(conda info --base 2>/dev/null || true)"
    fi
    celatlas_env_try_conda_root "$HOME/miniforge3"
    celatlas_env_try_conda_root "$HOME/mambaforge"
    celatlas_env_try_conda_root "$HOME/miniconda3"
    celatlas_env_try_conda_root "$HOME/anaconda3"
  fi

  if [[ -z "${CELATLAS_ENV_PREFIX:-}" && -n "$CELATLAS_CONDA_ROOT" ]]; then
    CELATLAS_ENV_PREFIX="${CELATLAS_CONDA_ROOT}/envs/${CELATLAS_ENV_NAME}"
  fi

  CELATLAS_AUTO_ACTIVATE="${CELATLAS_AUTO_ACTIVATE:-0}"
  CELATLAS_USE_ENV_PATH="${CELATLAS_USE_ENV_PATH:-0}"
}

celatlas_env_is_help_request() {
  for arg in "$@"; do
    case "$arg" in
      --help|-h)
        return 0
        ;;
    esac
  done
  return 1
}

celatlas_env_activate() {
  if [[ -z "${CELATLAS_CONDA_ROOT:-}" ]]; then
    echo "ERROR: CELATLAS_CONDA_ROOT is not set and conda root was not found"
    exit 1
  fi

  local conda_sh="${CELATLAS_CONDA_ROOT}/etc/profile.d/conda.sh"
  if [[ ! -f "$conda_sh" ]]; then
    echo "ERROR: conda.sh not found: ${conda_sh}"
    exit 1
  fi

  source "$conda_sh"
  conda activate "$CELATLAS_ENV_NAME" || {
    echo "ERROR: Failed to activate conda env '${CELATLAS_ENV_NAME}'"
    exit 1
  }

  if [[ -z "${CELATLAS_ENV_PREFIX:-}" ]]; then
    CELATLAS_ENV_PREFIX="${CONDA_PREFIX:-}"
  fi
}

celatlas_env_export_runtime() {
  if [[ -n "${CELATLAS_ENV_PREFIX:-}" && -x "${CELATLAS_ENV_PREFIX}/bin/python" ]]; then
    export STARDIST_PYTHON="${STARDIST_PYTHON:-${CELATLAS_ENV_PREFIX}/bin/python}"
  fi

  if [[ -n "${CELATLAS_ENV_PREFIX:-}" && -d "${CELATLAS_ENV_PREFIX}/lib" ]]; then
    case ":${LD_LIBRARY_PATH:-}:" in
      *":${CELATLAS_ENV_PREFIX}/lib:"*) ;;
      *) export LD_LIBRARY_PATH="${CELATLAS_ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}" ;;
    esac
  fi

  if [[ "$CELATLAS_USE_ENV_PATH" == "1" && -n "${CELATLAS_ENV_PREFIX:-}" && -d "${CELATLAS_ENV_PREFIX}/bin" ]]; then
    case ":${PATH:-}:" in
      *":${CELATLAS_ENV_PREFIX}/bin:"*) ;;
      *) export PATH="${CELATLAS_ENV_PREFIX}/bin:${PATH:-}" ;;
    esac
  fi

  export PYTHONNOUSERSITE="${PYTHONNOUSERSITE:-1}"
  export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR:-/tmp}/celatlas_mplconfig}"
  export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${TMPDIR:-/tmp}/celatlas_numba_cache}"

  if [[ "${CELATLAS_SKIP_ENV_MKDIR:-0}" != "1" ]]; then
    mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
  fi
}

celatlas_env_setup() {
  celatlas_env_init_defaults

  if [[ "$CELATLAS_AUTO_ACTIVATE" == "1" ]] && ! celatlas_env_is_help_request "$@"; then
    if [[ "${CONDA_DEFAULT_ENV:-}" != "$CELATLAS_ENV_NAME" ]]; then
      celatlas_env_activate
    elif [[ -z "${CELATLAS_ENV_PREFIX:-}" ]]; then
      CELATLAS_ENV_PREFIX="${CONDA_PREFIX:-}"
    fi
  fi

  if [[ -z "${CELATLAS_ENV_PREFIX:-}" && -n "${CONDA_PREFIX:-}" && "${CONDA_DEFAULT_ENV:-}" == "$CELATLAS_ENV_NAME" ]]; then
    CELATLAS_ENV_PREFIX="$CONDA_PREFIX"
  fi

  celatlas_env_export_runtime
}

celatlas_env_summary() {
  cat <<EOF
CELATLAS_REPO_DIR=${CELATLAS_REPO_DIR:-}
CELATLAS_CONFIG_FILE=${CELATLAS_CONFIG_FILE:-}
CELATLAS_ENV_NAME=${CELATLAS_ENV_NAME:-}
CELATLAS_ENV_PREFIX=${CELATLAS_ENV_PREFIX:-}
CELATLAS_AUTO_ACTIVATE=${CELATLAS_AUTO_ACTIVATE:-0}
CELATLAS_USE_ENV_PATH=${CELATLAS_USE_ENV_PATH:-0}
STARDIST_PYTHON=${STARDIST_PYTHON:-}
PYTHONNOUSERSITE=${PYTHONNOUSERSITE:-}
MPLCONFIGDIR=${MPLCONFIGDIR:-}
NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-}
EOF
}
