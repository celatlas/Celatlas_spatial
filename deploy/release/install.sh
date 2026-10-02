#!/usr/bin/env bash
set -euo pipefail

bundle_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
env_name="celatlas18"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda was not found. Install Miniforge or make conda available on PATH." >&2
  exit 1
fi

if command -v mamba >/dev/null 2>&1; then
  conda_tool=mamba
else
  conda_tool=conda
fi

if conda env list | awk '{print $1}' | grep -Fxq "$env_name"; then
  "$conda_tool" env update -n "$env_name" -f "$bundle_dir/envs/celatlas18.yml" --prune
else
  "$conda_tool" env create -f "$bundle_dir/envs/celatlas18.yml"
fi

"$conda_tool" install -y -n "$env_name" -c conda-forge -c bioconda subread=2.0.6 samtools=1.20
wheel_path="$(find "$bundle_dir/dist" -maxdepth 1 -name 'celatlas_spatial-*.whl' -print -quit)"
if [[ -z "$wheel_path" ]]; then
  echo "ERROR: Celatlas wheel is missing from $bundle_dir/dist" >&2
  exit 1
fi
conda run -n "$env_name" python -m pip install "$wheel_path" --no-deps

if [[ -d "$bundle_dir/models/stardist" ]]; then
  mkdir -p "$HOME/.keras/models/StarDist2D"
  cp -a "$bundle_dir/models/stardist/." "$HOME/.keras/models/StarDist2D/"
fi

if [[ ! -f "$bundle_dir/models/swin_tiny.pth" ]]; then
  echo "WARNING: swin_tiny.pth is not included. Copy it to the configured paths.src_dir before running binSegment." >&2
fi
echo "Installed Celatlas Spatial into conda environment: $env_name"
