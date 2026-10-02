#!/bin/bash
set -eE

# =========================
# Celatlas v1.8 Reanalysis Entrypoint
# =========================
# Purpose:
#   Re-run downstream analysis from existing 05.count outputs.
#   In HE workflows this is the preferred entrypoint after manual mask/
#   registration has been finalized.

# =========================
# OpenBLAS Thread Limit
# =========================
export OPENBLAS_NUM_THREADS=32
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export VECLIB_MAXIMUM_THREADS=32
export NUMEXPR_NUM_THREADS=32

# =========================
# Helpers
# =========================
run_command() {
  CELATLAS_LAST_COMMAND="$(printf '%q ' "$@")"
  echo "Running:" "$@"
  "$@"
  if [ $? -ne 0 ]; then
    echo "ERROR: Command failed: $*"
    exit 1
  fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/celatlas_env.sh"
celatlas_env_setup "$@"
source "${SCRIPT_DIR}/scripts/celatlas_runtime_paths.sh"
source "${SCRIPT_DIR}/scripts/celatlas_post_binsegment.sh"
source "${SCRIPT_DIR}/scripts/celatlas_failure_report.sh"
celatlas_post_binsegment_init_defaults
CELATLAS_WORKFLOW="reanalysis"

show_help() {
  cat <<'EOF'
=====================================================================
Celatlas v1.8 Reanalysis Script
=====================================================================

USAGE:
  ./Celatlas_reanalysis.sh <sample_name> <chip_number> <casno> <chemistry> <Species> <method> <mode> [options]
  ./Celatlas_reanalysis.sh <chip_number> <casno> <chemistry> <Species> <method> <mode> [options]

  ./Celatlas_reanalysis.sh --chip <chip> --casno <casno> --chemistry <chem> \
    --species <species> --method <method> --mode <mode> [options]

PARAMETERS:
  --sample, -s <name>        Sample name label, optional
  --chip, -c <chip>          Chip/slide number, e.g. SX000151_B1
  --casno <case>             Case/project number
  --chemistry <version>      Chemistry version, e.g. BBV2.4
  --species <species>        Species name, e.g. Mus_musculus or Homo_sapiens
  --method, -m <method>      image, ssDNA, gene_expr, or HE
  --mode <mode>              strna or scrna

REANALYSIS OPTIONS:
  --thread, -t <num>         CPU threads, default: auto-detected
  --bin, -b <sizes>          Analysis bin size(s), default: 50
  --pixelSize <um>           Pixel size in micrometers, default: 0.5
  --cell_num <num>           Expected cell number, default: 50000
  --no-clean                 Do not remove previous 06.segment/07.outs
  --skip-binsegment          Reuse existing 06.segment/01.binsegment
  --skip-analysis            Skip 07.outs analysis
  --skip-report              Skip report generation
  --gene-mask-filter         Restrict gene_expr or ssDNA/image binning to the
                             GEM/manual gene mask region

PATH OPTIONS:
  --workspace <path>         Workspace root, default: $CELATLAS_WORKSPACE
  --reference_dir <path>     Reference root, default: <workspace>/reference
  --genomeDir <path>         Explicit genome directory
  --mask_dir <path>          Source barcode/mask directory, default: <workspace>/ST_mask
  --image_dir <path>         Source image directory, default: <workspace>/images
  --sampledir <path>         Result directory, default: <workspace>/results/<casno>/<chip>
  --src_dir <path>           Source/model directory, default: <workspace>/src

POST-MANUAL-REGISTRATION MODULES:
EOF
  celatlas_post_binsegment_help
  cat <<'EOF'
EXAMPLES:
  # HE reanalysis after manual registration, plus cavity + StarDist
  ./Celatlas_reanalysis.sh SX000151_B1 SX_test BBV2.4 Mus_musculus HE strna \
    --thread 32 \
    --bin 50 \
    --enable-cavity-filter \
    --enable-stardist-cell-segment

  # Only run post-binSegment modules on an existing final 06.segment/01.binsegment
  ./Celatlas_reanalysis.sh SX000151_B1 SX_test BBV2.4 Mus_musculus HE strna \
    --skip-binsegment --skip-analysis --skip-report \
    --enable-cavity-filter --enable-stardist-cell-segment

EOF
}

# =========================
# Input Parsing
# =========================
sample_name=""
chip_number=""
casno=""
chemistry=""
Species=""
method=""
mode=""
workspace_dir="${CELATLAS_WORKSPACE}"
reference_dir=""
genomeDir=""
sampledir=""
src_dir=""
mask_dir=""
image_dir=""
OPT_THREAD=""
OPT_BIN=""
OPT_PIXELSIZE=""
OPT_CELL_NUM=""
OPT_NO_CLEAN=0
OPT_SKIP_BINSEGMENT=0
OPT_SKIP_ANALYSIS=0
OPT_SKIP_REPORT=0
OPT_GENE_MASK_FILTER=0

positional_args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      show_help
      exit 0
      ;;
    --sample|-s)
      sample_name="$2"
      shift 2
      ;;
    --chip|-c)
      chip_number="$2"
      shift 2
      ;;
    --casno)
      casno="$2"
      shift 2
      ;;
    --chemistry)
      chemistry="$2"
      shift 2
      ;;
    --species)
      Species="$2"
      shift 2
      ;;
    --method|-m)
      method="$2"
      shift 2
      ;;
    --mode)
      mode="$2"
      shift 2
      ;;
    --workspace)
      workspace_dir="$2"
      shift 2
      ;;
    --reference_dir)
      reference_dir="$2"
      shift 2
      ;;
    --genomeDir)
      genomeDir="$2"
      shift 2
      ;;
    --mask_dir)
      mask_dir="$2"
      shift 2
      ;;
    --image_dir)
      image_dir="$2"
      shift 2
      ;;
    --sampledir)
      sampledir="$2"
      shift 2
      ;;
    --src_dir)
      src_dir="$2"
      shift 2
      ;;
    --thread|-t)
      OPT_THREAD="$2"
      shift 2
      ;;
    --bin|-b)
      OPT_BIN="$2"
      shift 2
      ;;
    --pixelSize)
      OPT_PIXELSIZE="$2"
      shift 2
      ;;
    --cell_num)
      OPT_CELL_NUM="$2"
      shift 2
      ;;
    --no-clean)
      OPT_NO_CLEAN=1
      shift
      ;;
    --skip-binsegment)
      OPT_SKIP_BINSEGMENT=1
      shift
      ;;
    --skip-analysis)
      OPT_SKIP_ANALYSIS=1
      shift
      ;;
    --skip-report)
      OPT_SKIP_REPORT=1
      shift
      ;;
    --gene-mask-filter)
      OPT_GENE_MASK_FILTER=1
      shift
      ;;
    --enable-cavity-filter|--disable-cavity-filter|--skip-cavity-filter|\
    --cavity-bins|--cavity-apply|--cavity-dry-run|--cavity-skip-qc-images|--cavity-mask-preset|\
    --enable-stardist-cell-segment|--enable-stardist|\
    --disable-stardist-cell-segment|--skip-stardist|--disable-stardist|\
    --stardist-python|--stardist-labels|--stardist-image|--stardist-tissue-bbox|--stardist-model|--stardist-model-dir|\
    --stardist-prob-thresh|--stardist-nms-thresh|--stardist-max-dim|\
    --stardist-scale|--stardist-n-tiles|--stardist-expand-pixels|\
    --stardist-min-umi|--stardist-min-genes|--stardist-skip-counts|--stardist-no-label-output|\
    --stardist-tiled-inference|--stardist-tile-size|--stardist-tile-overlap|\
    --stardist-tile-merge-overlap|--stardist-fluorescence-channel|\
    --stardist-tile-cache-dir|--stardist-exclude-rectangle)
      celatlas_post_binsegment_parse_option "$@"
      shift "$CELATLAS_POST_ARG_SHIFT"
      ;;
    -*)
      echo "ERROR: Unknown option '$1'"
      echo "Run '$0 --help' for usage information"
      exit 1
      ;;
    *)
      positional_args+=("$1")
      shift
      ;;
  esac
done

if [ ${#positional_args[@]} -eq 7 ]; then
  sample_name="${positional_args[0]}"
  chip_number="${positional_args[1]}"
  casno="${positional_args[2]}"
  chemistry="${positional_args[3]}"
  Species="${positional_args[4]}"
  method="${positional_args[5]}"
  mode="${positional_args[6]}"
elif [ ${#positional_args[@]} -eq 6 ]; then
  chip_number="${positional_args[0]}"
  casno="${positional_args[1]}"
  chemistry="${positional_args[2]}"
  Species="${positional_args[3]}"
  method="${positional_args[4]}"
  mode="${positional_args[5]}"
elif [ ${#positional_args[@]} -gt 0 ]; then
  echo "ERROR: Invalid number of positional arguments (${#positional_args[@]})"
  echo "Expected 6 or 7 positional arguments, or use named arguments."
  exit 1
fi

if [[ -z "$chip_number" || -z "$casno" || -z "$chemistry" || -z "$Species" || -z "$method" || -z "$mode" ]]; then
  echo "ERROR: Missing required parameters"
  echo "Run '$0 --help' for usage information"
  exit 1
fi

sample="$chip_number"
export CELATLAS_SAMPLE_NAME="$sample_name"
export CELATLAS_CHIP_NUMBER="$chip_number"

if [[ "$method" != "image" && "$method" != "ssDNA" && "$method" != "gene_expr" && "$method" != "HE" ]]; then
  echo "ERROR: Invalid method '$method'. Must be image, ssDNA, gene_expr, or HE."
  exit 1
fi
if [[ "$method" == "ssDNA" ]]; then
  echo "Note: Using ssDNA mode (equivalent to image mode)"
  method="image"
fi
celatlas_post_binsegment_finalize_defaults "$method"
if [[ "$mode" != "scrna" && "$mode" != "strna" ]]; then
  echo "ERROR: Invalid mode '$mode'. Must be scrna or strna."
  exit 1
fi

# =========================
# Runtime Configuration
# =========================
if command -v nproc >/dev/null 2>&1; then
  CPU_CORES=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
  CPU_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo 16)
else
  CPU_CORES=16
fi
thread=${OPT_THREAD:-${THREADS:-$CPU_CORES}}

if [ "$thread" -gt 32 ]; then
  blas_threads=32
  echo "Note: Limiting BLAS threads to 32 (user thread count: $thread)"
else
  blas_threads=$thread
fi

export OMP_NUM_THREADS=$thread
export OPENBLAS_NUM_THREADS=$blas_threads
export MKL_NUM_THREADS=$blas_threads
export VECLIB_MAXIMUM_THREADS=$blas_threads
export BLIS_NUM_THREADS=$blas_threads
export NUMEXPR_MAX_THREADS=$thread
export NUMEXPR_NUM_THREADS=$thread
export BINSEGMENT_MAX_PARALLEL=${BINSEGMENT_MAX_PARALLEL:-2}

bin=${OPT_BIN:-50}
pixelSize=${OPT_PIXELSIZE:-0.5}
cell_num=${OPT_CELL_NUM:-50000}
feature_type="gene"
reference_dir=${reference_dir:-"${CELATLAS_REFERENCE_DIR}"}
src_dir=${src_dir:-"${CELATLAS_SRC_DIR}"}
sampledir=${sampledir:-"${CELATLAS_RESULTS_ROOT}/${casno}/${sample}"}
mask_dir=${mask_dir:-"${CELATLAS_MASK_DIR}"}
image_dir=${image_dir:-"${CELATLAS_IMAGE_DIR}"}
runtime_mask_dir=""
if [[ "$mode" != "scrna" ]]; then
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "${sampledir}")"
fi
genomeDir=${genomeDir:-"${reference_dir}/${Species}"}

if [[ "$mode" == "scrna" ]]; then
  mkdir -p "${sampledir}"
else
  mkdir -p "${runtime_mask_dir}" "${sampledir}"
fi
log_file="${sampledir}/reanalysis_pipeline.log"
exec > >(tee -a "$log_file") 2>&1
trap 'celatlas_on_pipeline_error "$?" "$LINENO" "$BASH_COMMAND"' ERR
trap 'celatlas_on_pipeline_exit "$?" "$LINENO" "$BASH_COMMAND"' EXIT

echo "======================================"
echo "Starting Celatlas v1.8 reanalysis at $(date)"
echo "======================================"
echo "sample=${sample} casno=${casno} chemistry=${chemistry} Species=${Species} method=${method} mode=${mode}"
echo "workspace_dir=${workspace_dir}"
echo "sampledir=${sampledir}"
echo "runtime_mask_dir=${runtime_mask_dir:-N/A}"
echo "mask_dir=${mask_dir}"
echo "image_dir=${image_dir}"
echo "genomeDir=${genomeDir}"
echo "threads=${thread} blas_threads=${blas_threads} bin=${bin} pixelSize=${pixelSize}"
echo "skip_binsegment=${OPT_SKIP_BINSEGMENT} skip_analysis=${OPT_SKIP_ANALYSIS} skip_report=${OPT_SKIP_REPORT} no_clean=${OPT_NO_CLEAN}"
echo "enable_cavity_filter=${OPT_ENABLE_CAVITY_FILTER} enable_stardist_cell_segment=${OPT_ENABLE_STARDIST_CELL_SEGMENT}"
echo "gene_mask_filter=${OPT_GENE_MASK_FILTER}"

# =========================
# Prerequisite Checks
# =========================
echo ""
echo "[Prerequisite Check]"

if [[ ! -d "${sampledir}/05.count" ]]; then
  echo "ERROR: 05.count directory not found: ${sampledir}/05.count"
  exit 1
fi
echo "OK: Found 05.count"

if [[ "$mode" != "scrna" && ! -f "${sampledir}/05.count/${sample}_count_detail.txt" ]]; then
  echo "ERROR: count_detail not found: ${sampledir}/05.count/${sample}_count_detail.txt"
  exit 1
fi
if [[ "$mode" != "scrna" ]]; then
  echo "OK: Found count_detail ($(du -h "${sampledir}/05.count/${sample}_count_detail.txt" | cut -f1))"
fi

if [[ "$mode" != "scrna" ]]; then
  if [[ ! -f "${runtime_mask_dir}/${sample}_FilterBarcodes.csv" ]]; then
    echo "[Setup] Runtime mask/input files missing; staging from source directories."
    celatlas_stage_spatial_inputs "$sample" "$sampledir" "$mask_dir" "$image_dir" "$method"
  else
    celatlas_stage_gene_masks "$sample" "$sampledir" "$mask_dir"
  fi

  required_files=(
    "${runtime_mask_dir}/${sample}_FilterBarcodes.csv"
  )
  if [[ "$method" == "image" ]]; then
    required_files+=("${runtime_mask_dir}/${sample}.tif")
  fi
  for f in "${required_files[@]}"; do
    if [[ ! -f "$f" ]]; then
      echo "ERROR: Required runtime input not found: $f"
      exit 1
    fi
    echo "OK: Found $(basename "$f")"
  done
  if [[ "$method" == "HE" ]]; then
    if ! celatlas_find_runtime_he_image "$sample" "$sampledir" "$image_dir" >/dev/null; then
      echo "ERROR: HE mode requires ${runtime_mask_dir}/${sample}_he.(tif|tiff|png|jpg|jpeg)"
      exit 1
    fi
    echo "OK: Found HE image"
  fi
fi

if [[ ! -d "${genomeDir}" ]]; then
  echo "ERROR: Reference genome not found: ${genomeDir}"
  exit 1
fi
echo "OK: Found reference genome"

ulimit -n 10240

# =========================
# Cleanup
# =========================
if [[ "$OPT_NO_CLEAN" != "1" ]]; then
  echo ""
  echo "[Cleanup]"
  if [[ "$OPT_SKIP_BINSEGMENT" != "1" && -d "${sampledir}/06.segment/01.binsegment" ]]; then
    echo "Removing ${sampledir}/06.segment/01.binsegment/"
    rm -rf "${sampledir}/06.segment/01.binsegment/"
  fi
  if [[ ( "$OPT_SKIP_BINSEGMENT" != "1" || "$OPT_ENABLE_STARDIST_CELL_SEGMENT" == "1" ) && -d "${sampledir}/06.segment/02.cellsegment" ]]; then
    echo "Removing ${sampledir}/06.segment/02.cellsegment/"
    rm -rf "${sampledir}/06.segment/02.cellsegment/"
  fi
  if [[ "$OPT_SKIP_ANALYSIS" != "1" && -d "${sampledir}/07.outs" ]]; then
    echo "Removing ${sampledir}/07.outs/"
    rm -rf "${sampledir}/07.outs/"
  fi
  if [[ "$OPT_SKIP_REPORT" != "1" && -f "${sampledir}/${sample}_spatial_analysis_report.html" ]]; then
    echo "Removing old report: ${sample}_spatial_analysis_report.html"
    rm -f "${sampledir}/${sample}_spatial_analysis_report.html"
  fi
else
  echo ""
  echo "[Cleanup] Skipped by --no-clean"
fi

# =========================
# binSegment Branches
# =========================
run_image_branch() {
  echo "======================================"
  echo "[Image Mode] Running 06.segment/01.binsegment"
  echo "======================================"
  local gene_mask_args=()
  if [[ "$OPT_GENE_MASK_FILTER" == "1" ]]; then
    gene_mask_args+=(--gene-mask-filter)
  fi
  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.segment/01.binsegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${runtime_mask_dir}" \
    --segment --tif "${runtime_mask_dir}/${sample}.tif" --model "${src_dir}/swin_tiny.pth" --method "image" \
    --count --count_detail "${sampledir}/05.count/${sample}_count_detail.txt" \
    "${gene_mask_args[@]}"
}

run_gene_expr_branch() {
  echo "======================================"
  echo "[Gene Expression Mode] Running 06.segment/01.binsegment"
  echo "======================================"
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'

  echo "gem_bin_size=${gem_bin_size} umi_min_threshold=${umi_threshold}"
  local gene_mask_args=()
  if [[ "$OPT_GENE_MASK_FILTER" == "1" ]]; then
    gene_mask_args+=(--gene-mask-filter)
  fi
  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.segment/01.binsegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${runtime_mask_dir}" \
    --segment --model "${src_dir}/swin_tiny.pth" --method "gene_expr" --count \
    --count_detail "${sampledir}/05.count/${sample}_count_detail.txt" \
    --gem-bin-size "${gem_bin_size}" \
    --umi-min-threshold "${umi_threshold}" \
    --enhance-params "${enhance_params}" \
    "${gene_mask_args[@]}"
}

run_he_branch() {
  echo "======================================"
  echo "[HE Mode] Running final 06.segment/01.binsegment after manual registration/mask update"
  echo "======================================"
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'
  local registration_type="${REGISTRATION_TYPE:-affine}"
  local he_image_found=""

  he_image_found="$(celatlas_find_runtime_he_image "$sample" "$sampledir" "$image_dir" || true)"
  if [[ -z "$he_image_found" ]]; then
    echo "ERROR: HE image not found under ${runtime_mask_dir}"
    exit 1
  fi

  echo "gem_bin_size=${gem_bin_size} umi_min_threshold=${umi_threshold} registration_type=${registration_type}"
  echo "HE image: ${he_image_found}"
  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.segment/01.binsegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${runtime_mask_dir}" \
    --segment --model "${src_dir}/swin_tiny.pth" --method "HE" --count \
    --count_detail "${sampledir}/05.count/${sample}_count_detail.txt" \
    --tif "${he_image_found}" \
    --gem-bin-size "${gem_bin_size}" \
    --umi-min-threshold "${umi_threshold}" \
    --enhance-params "${enhance_params}" \
    --registration-type "${registration_type}"
}

run_scrna_reanalysis() {
  echo "======================================"
  echo "scRNA Reanalysis"
  echo "======================================"
  matrix_dir="${sampledir}/05.count/${sample}_filtered_feature_bc_matrix"
  if [[ ! -d "$matrix_dir" ]]; then
    echo "ERROR: Filtered matrix not found: ${matrix_dir}"
    exit 1
  fi
  if [[ "$OPT_NO_CLEAN" != "1" && -d "${sampledir}/06_analysis_wrapper" ]]; then
    rm -rf "${sampledir}/06_analysis_wrapper"
  fi
  if [[ "$OPT_SKIP_ANALYSIS" != "1" ]]; then
    run_command celatlas_spatial rna analysis \
      --outdir "${sampledir}/06_analysis_wrapper" \
      --sample "${sample}" \
      --thread "${thread}" \
      --genomeDir "${genomeDir}" \
      --matrix_file "${matrix_dir}" \
      --assay scrna
  fi
}

run_spatial_analysis() {
  if [[ "$OPT_SKIP_ANALYSIS" == "1" ]]; then
    echo "Step 07.outs skipped by --skip-analysis"
    return 0
  fi

  echo "======================================"
  echo "Step 07: Outs Analysis"
  echo "======================================"

  local skip_tsne_flag=""
  if [[ "$bin" == "10" || "$bin" == "20" ]]; then
    skip_tsne_flag="--skip-tsne"
    echo "[Dimensional Reduction] Bin size ${bin}: using UMAP only"
  else
    echo "[Dimensional Reduction] Bin size ${bin}: using t-SNE and UMAP"
  fi

  run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.outs/binned_outputs" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --square_bin_dir "${sampledir}/06.segment/01.binsegment/square_bin" \
    --pixel-size "${pixelSize}" --bin "${bin}" ${skip_tsne_flag}

  local cell_matrix_dir="${sampledir}/06.segment/02.cellsegment/cell_matrix"
  if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" == "1" && -d "${cell_matrix_dir}" && -f "${cell_matrix_dir}/matrix.mtx.gz" ]]; then
    echo "======================================"
    echo "Step 07: Cell-Segmented Analysis"
    echo "======================================"
    run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.outs/cellsegmented_outputs" --sample "${sample}" \
      --thread "${thread}" --genomeDir "${genomeDir}" --assay cell \
      --cell-seg-dir "${sampledir}/06.segment/02.cellsegment" \
      --skip-tsne
  fi
}

run_spatial_report() {
  if [[ "$OPT_SKIP_REPORT" == "1" ]]; then
    echo "Spatial report skipped by --skip-report"
    return 0
  fi

  echo "======================================"
  echo "Step 08: Spatial Report"
  echo "======================================"
  local out_html="${sample}_spatial_analysis_report.html"
  local cell_report_args=()
  if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" != "1" ]]; then
    cell_report_args+=(--exclude-cell-segmentation)
  fi
  if command -v celatlas_spatial_report >/dev/null 2>&1; then
    run_command celatlas_spatial_report \
      "${sampledir}" "${sample}" \
      --chemistry "${chemistry}" \
      --species "${Species}" \
      --method "${method}" \
      --output-filename "${out_html}" \
      --verbose \
      "${cell_report_args[@]}"
  else
    local report_generator="${CELATLAS_SPATIAL_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/spatial_report_generator.py}"
    if [[ -f "$report_generator" ]]; then
      run_command python3 "$report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --method "${method}" \
        --output-filename "${out_html}" \
        --verbose \
        "${cell_report_args[@]}"
    else
      echo "WARNING: Spatial report generator not found. Skipping report."
    fi
  fi
}

run_scrna_report() {
  if [[ "$OPT_SKIP_REPORT" == "1" ]]; then
    echo "scRNA report skipped by --skip-report"
    return 0
  fi

  local out_html="${sample}_scrna_analysis_report.html"
  if command -v celatlas_scrna_report >/dev/null 2>&1; then
    run_command celatlas_scrna_report \
      "${sampledir}" "${sample}" \
      --chemistry "${chemistry}" \
      --species "${Species}" \
      --output-filename "${out_html}" \
      --verbose
  else
    local report_generator="${CELATLAS_SCRNA_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/scrna_report_generator.py}"
    if [[ -f "$report_generator" ]]; then
      run_command python3 "$report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --output-filename "${out_html}" \
        --verbose
    else
      echo "WARNING: scRNA report generator not found. Skipping report."
    fi
  fi
}

# =========================
# Main
# =========================
if [[ "$mode" == "scrna" ]]; then
  run_scrna_reanalysis
  run_scrna_report
else
  if [[ "$OPT_SKIP_BINSEGMENT" == "1" ]]; then
    echo "Step 06.segment/01.binsegment skipped by --skip-binsegment"
    if [[ ! -d "${sampledir}/06.segment/01.binsegment" ]]; then
      echo "ERROR: Existing 06.segment/01.binsegment not found: ${sampledir}/06.segment/01.binsegment"
      exit 1
    fi
  else
    if [[ "$method" == "image" ]]; then
      run_image_branch
    elif [[ "$method" == "HE" ]]; then
      run_he_branch
    else
      run_gene_expr_branch
    fi
    celatlas_sync_gem_expression_to_mask_dir "$sampledir"
  fi

  celatlas_run_post_binsegment_optional_modules
  run_spatial_analysis
  run_spatial_report
fi

echo "======================================"
echo "Celatlas v1.8 reanalysis completed at $(date)"
echo "======================================"
