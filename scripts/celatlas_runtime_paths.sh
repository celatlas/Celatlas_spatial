#!/usr/bin/env bash

# Shared runtime path helpers for Celatlas shell entrypoints.
# Source this file; do not execute it directly.

celatlas_mask_dir_for_sample() {
  local sampledir="$1"
  printf '%s/06.segment/mask' "$sampledir"
}

celatlas_stage_manual_masks() {
  local sample="$1"
  local sampledir="$2"
  local source_mask_dir="$3"

  local runtime_mask_dir
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "$sampledir")"
  mkdir -p "$runtime_mask_dir"

  local mask_name
  for mask_name in \
    "manual_mask.png" \
    "${sample}_manual_mask.png" \
    "${sample}_mask_manual.png" \
    "mask_manual.png" \
    "he_manual_mask.png" \
    "${sample}_he_manual_mask.png" \
    "${sample}_he_mask_manual.png" \
    "he_mask_manual.png"; do
    if [[ -f "${source_mask_dir}/${mask_name}" ]]; then
      cp -n "${source_mask_dir}/${mask_name}" "$runtime_mask_dir/"
      echo "[Setup] Copied manual mask: ${mask_name}"
    fi
  done
}

celatlas_stage_spatial_inputs() {
  local sample="$1"
  local sampledir="$2"
  local source_mask_dir="$3"
  local source_image_dir="$4"
  local method="$5"

  local runtime_mask_dir
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "$sampledir")"
  mkdir -p "$runtime_mask_dir"

  local filter_barcodes_csv="${source_mask_dir}/${sample}_FilterBarcodes.csv"
  local tissue_bbox_csv="${source_mask_dir}/${sample}_tissue_bbox.csv"
  local barcode_h5="${source_mask_dir}/${sample}.barcodeToPos.h5"

  if [[ ! -f "$filter_barcodes_csv" ]]; then
    echo ""
    echo "ERROR: spatial barcode positions file not found."
    echo "       Required:"
    echo "       - $filter_barcodes_csv"
    echo "       Optional legacy tissue bbox CSV:"
    echo "       - $tissue_bbox_csv"
    echo "       Optional barcode whitelist/canvas H5:"
    echo "       - $barcode_h5"
    echo ""
    return 1
  fi

  cp -n "$filter_barcodes_csv" "$runtime_mask_dir/"
  if [[ -f "$tissue_bbox_csv" ]]; then
    cp -n "$tissue_bbox_csv" "$runtime_mask_dir/"
  fi
  if [[ -f "$barcode_h5" ]]; then
    cp -n "$barcode_h5" "$runtime_mask_dir/"
  fi

  celatlas_stage_manual_masks "$sample" "$sampledir" "$source_mask_dir"

  local he_found=""
  for ext in tif tiff png jpg jpeg; do
    if [[ -f "${source_image_dir}/${sample}_he.${ext}" ]]; then
      he_found="${source_image_dir}/${sample}_he.${ext}"
      cp -n "$he_found" "$runtime_mask_dir/"
      break
    fi
  done

  if [[ "$method" == "HE" && -z "$he_found" ]]; then
    echo ""
    echo "ERROR: HE mode requires H&E staining image."
    echo "       Missing: ${source_image_dir}/${sample}_he.(tif|tiff|png|jpg|jpeg)"
    echo ""
    return 1
  fi

  if [[ "$method" == "image" ]]; then
    local tissue_image="${source_image_dir}/${sample}.tif"
    if [[ ! -f "$tissue_image" ]]; then
      echo ""
      echo "ERROR: image/ssDNA mode requires tissue image file."
      echo "       Missing: $tissue_image"
      echo ""
      return 1
    fi
    cp -n "$tissue_image" "$runtime_mask_dir/"
  fi

  echo "[Setup] Runtime mask/input directory: $runtime_mask_dir"
}

celatlas_find_runtime_he_image() {
  local sample="$1"
  local sampledir="$2"
  local image_dir="${3:-}"

  local runtime_mask_dir
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "$sampledir")"
  local ext path
  for ext in tif tiff png jpg jpeg; do
    path="${runtime_mask_dir}/${sample}_he.${ext}"
    if [[ -f "$path" ]]; then
      printf '%s\n' "$path"
      return 0
    fi
    if [[ -n "$image_dir" ]]; then
      path="${image_dir}/${sample}_he.${ext}"
      if [[ -f "$path" ]]; then
        printf '%s\n' "$path"
        return 0
      fi
    fi
  done
  return 1
}

celatlas_sync_gem_expression_to_mask_dir() {
  local sampledir="$1"
  local source="${sampledir}/06.segment/01.binsegment/images/overlays/1_gem_expression.png"
  local runtime_mask_dir
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "$sampledir")"
  mkdir -p "$runtime_mask_dir"

  if [[ -f "$source" ]]; then
    cp -f "$source" "${runtime_mask_dir}/1_gem_expression.png"
    echo "[Setup] Copied GEM expression overlay to ${runtime_mask_dir}/1_gem_expression.png"
  else
    echo "WARNING: GEM expression overlay not found for mask editing: $source"
  fi
}

celatlas_fastq_manifest_status() {
  local previous="$1"
  local current="$2"

  if [[ ! -f "$previous" ]]; then
    printf 'new\n'
    return 0
  fi
  if cmp -s "$current" "$previous"; then
    printf 'same\n'
    return 0
  fi
  if awk 'NR==FNR { current[$0]=1; next } !($0 in current) { missing=1 } END { exit missing ? 1 : 0 }' "$current" "$previous"; then
    local extra_count
    extra_count="$(awk 'NR==FNR { previous[$0]=1; next } !($0 in previous) { count++ } END { print count + 0 }' "$previous" "$current")"
    if (( extra_count > 0 )); then
      printf 'added\n'
    else
      printf 'reordered\n'
    fi
    return 0
  fi
  printf 'changed\n'
}

celatlas_has_fastq_dependent_outputs() {
  local sampledir="$1"
  local sample="$2"
  local path
  for path in \
    "${sampledir}/00.sample" \
    "${sampledir}/01.barcode" \
    "${sampledir}/02.cutadapt" \
    "${sampledir}/03.star" \
    "${sampledir}/04.featureCounts" \
    "${sampledir}/05.count" \
    "${sampledir}/06.segment/01.binsegment" \
    "${sampledir}/06.segment/02.cellsegment" \
    "${sampledir}/07.outs" \
    "${sampledir}/06_analysis_wrapper" \
    "${sampledir}/${sample}_spatial_analysis_report.html" \
    "${sampledir}/${sample}_scrna_analysis_report.html" \
    "${sampledir}/${sample}_report.html"; do
    if [[ -e "$path" ]]; then
      return 0
    fi
  done
  return 1
}

celatlas_reset_one_path() {
  local archive_root="$1"
  local path="$2"

  if [[ ! -e "$path" ]]; then
    return 1
  fi

  if [[ "${CELATLAS_RERUN_ARCHIVE:-0}" == "1" ]]; then
    mkdir -p "$archive_root"
    local base target index
    base="$(basename "$path")"
    target="${archive_root}/${base}"
    index=1
    while [[ -e "$target" ]]; do
      target="${archive_root}/${base}.${index}"
      index=$((index + 1))
    done
    mv "$path" "$target"
    echo "[FASTQ] Archived stale output: $path -> $target"
  else
    rm -rf "$path"
    echo "[FASTQ] Removed stale output: $path"
  fi
  return 0
}

celatlas_reset_fastq_dependent_outputs() {
  local sampledir="$1"
  local sample="$2"
  local reason="$3"
  local timestamp archive_root reset path

  timestamp="$(date +%Y%m%d_%H%M%S)_$$"
  archive_root="${sampledir}/.rerun_archive/fastq_${reason}_${timestamp}"
  reset=0

  for path in \
    "${sampledir}/00.sample" \
    "${sampledir}/01.barcode" \
    "${sampledir}/02.cutadapt" \
    "${sampledir}/03.star" \
    "${sampledir}/04.featureCounts" \
    "${sampledir}/05.count" \
    "${sampledir}/06.segment/01.binsegment" \
    "${sampledir}/06.segment/02.cellsegment" \
    "${sampledir}/07.outs" \
    "${sampledir}/06_analysis_wrapper" \
    "${sampledir}/${sample}_spatial_analysis_report.html" \
    "${sampledir}/${sample}_scrna_analysis_report.html" \
    "${sampledir}/${sample}_report.html"; do
    if celatlas_reset_one_path "$archive_root" "$path"; then
      reset=1
    fi
  done

  if (( reset == 1 )); then
    echo "[FASTQ] Preserved manual mask/runtime inputs under ${sampledir}/06.segment/mask"
  elif [[ "${CELATLAS_RERUN_ARCHIVE:-0}" == "1" ]]; then
    rmdir "$archive_root" 2>/dev/null || true
  fi
}
