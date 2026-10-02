#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
version="$(sed -n 's/^__VERSION__ = "\([^"]*\)"/\1/p' "$repo_dir/celatlas_spatial/__init__.py")"
if [[ -z "$version" ]]; then
  echo "ERROR: Could not determine Celatlas version." >&2
  exit 1
fi

output_dir="${1:-$(dirname "$repo_dir")/celatlas-spatial-v${version}-linux-x86_64}"
archive_path="${output_dir}.tar.gz"
if [[ -e "$output_dir" || -e "$archive_path" ]]; then
  echo "ERROR: Output already exists: $output_dir or $archive_path" >&2
  exit 1
fi

mkdir -p "$output_dir"/{dist,docs,envs,config-templates,models,release}
rm -rf "$repo_dir/build" "$repo_dir/celatlas_spatial.egg-info"
python -m build --no-isolation --wheel --sdist --outdir "$output_dir/dist" "$repo_dir"
python "$repo_dir/scripts/check_public_release.py" "$output_dir"/dist/celatlas_spatial-*.whl "$output_dir"/dist/celatlas_spatial-*.tar.gz

cp "$repo_dir/README.md" "$repo_dir/LICENSE.txt" "$output_dir/"
cp -a "$repo_dir/envs/." "$output_dir/envs/"
cp "$repo_dir/docs/celatlas_spatial_manual_zh.md" \
  "$repo_dir/docs/celatlas_spatial_manual_en.md" \
  "$repo_dir/docs/linux_deployment_guide.md" "$output_dir/docs/"
cp "$repo_dir/configs/runner.yaml.example" "$repo_dir/configs/celatlas.env.example" "$output_dir/config-templates/"
if [[ -d "$repo_dir/demo_data" ]]; then
  cp -a "$repo_dir/demo_data" "$output_dir/"
fi
cp "$repo_dir/deploy/release/install.sh" "$repo_dir/deploy/release/verify.sh" "$repo_dir/deploy/release/MISSING_RUNTIME_ASSETS.md" "$output_dir/release/"
cp "$repo_dir/scripts/create_fastq_demo.py" "$output_dir/release/create_fastq_demo.py"
chmod 0755 "$output_dir/release/create_fastq_demo.py"
chmod 0755 "$output_dir/release/install.sh" "$output_dir/release/verify.sh"

model_source="$repo_dir/.stardist_keras/models/StarDist2D"
if [[ -d "$model_source/2D_versatile_he/2D_versatile_he_extracted" ]]; then
  mkdir -p "$output_dir/models/stardist"
  cp -a "$model_source/2D_versatile_he" "$output_dir/models/stardist/"
fi

swin_model="${CELATLAS_SWIN_MODEL:-$repo_dir/swin_tiny.pth}"
if [[ -f "$swin_model" ]]; then
  cp "$swin_model" "$output_dir/models/swin_tiny.pth"
  swin_model_note="Included from ${swin_model}; SHA256 $(sha256sum "$swin_model" | awk '{print $1}')."
else
  swin_archive="$repo_dir/swin_tiny_model.tar.gz"
  if [[ -f "$swin_archive" ]]; then
    tar -xOf "$swin_archive" --wildcards '*/swin_tiny.pth' > "$output_dir/models/swin_tiny.pth"
    swin_model_note="Included Swin tissue-segmentation model; SHA256 $(sha256sum "$output_dir/models/swin_tiny.pth" | awk '{print $1}'). Copy it into the configured paths.src_dir before running binSegment."
  else
    swin_model_note="Not included; configure paths.src_dir with a swin_tiny.pth model before running binSegment."
  fi
fi

git_revision="$(git -C "$repo_dir" rev-parse --short HEAD 2>/dev/null || printf 'not-a-git-checkout')"
modified_count="$(git -C "$repo_dir" status --short 2>/dev/null | wc -l | tr -d ' ')"
cat > "$output_dir/RELEASE_NOTES.md" <<EOF
# Celatlas Spatial v${version} Linux Deployment Bundle

Built from repository revision: ${git_revision}.

This bundle was created from the current working tree, which had ${modified_count}
uncommitted path(s) when packaged. Treat it as a traceable deployment snapshot,
not a clean source-control release, until the source revision is committed and tagged.

Install with:

./release/install.sh
./release/verify.sh

See release/MISSING_RUNTIME_ASSETS.md before an analysis run.

Swin tissue-segmentation model: ${swin_model_note}
EOF

(cd "$output_dir" && find . -type f ! -name CHECKSUMS.sha256 -print0 | sort -z | xargs -0 sha256sum > CHECKSUMS.sha256)
tar -C "$(dirname "$output_dir")" -czf "$archive_path" "$(basename "$output_dir")"
printf 'Deployment folder: %s\nArchive: %s\n' "$output_dir" "$archive_path"
