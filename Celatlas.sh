#!/bin/bash
set -e

# =========================
# OpenBLAS Thread Limit
# =========================
# OpenBLAS Thread Limit (Initial)
# Initial conservative setting - will be adjusted based on CPU count later
# For high-resolution bins (10/20), t-SNE will be skipped automatically
export OPENBLAS_NUM_THREADS=16
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16

# =========================
# Helpers
# =========================
CLEAN_INTERMEDIATE=${CLEAN_INTERMEDIATE:-1}

safe_rm() {
  local f="$1"
  if [[ -f "$f" ]]; then rm -f "$f" && echo "[CLEAN] removed $f"; fi
}

safe_glob_rm() {
  local pattern="$1"
  shopt -s nullglob
  local files=( $pattern )
  if (( ${#files[@]} )); then
    rm -f "${files[@]}" && echo "[CLEAN] removed ${#files[@]} files in pattern $pattern"
  fi
  shopt -u nullglob
}

run_command () {
  echo "Running:" "$@"
  "$@"
}

# =========================
# Inputs - Hybrid Mode (Positional + Named Arguments)
# =========================
# Usage Examples:
#   Legacy mode (positional):
#     $0 <sample_name> <chip_number> <casno> <chemistry> <Species> <method> <mode>
#     $0 <chip_number> <casno> <chemistry> <Species> <method> <mode>
#
#   Modern mode (named arguments):
#     $0 --sample DEMO --chip ST110250_A1 --casno HE_test --chemistry BBV2.4 \
#        --species Mus_musculus --method HE --mode strna
#
#   Hybrid mode (mix both):
#     $0 DEMO ST110250_A1 HE_test BBV2.4 Mus_musculus HE strna --thread 32 --bin 10,50,100
#
#   Optional parameters (use defaults if not specified):
#     --thread <num>           Number of CPU threads (default: auto-detect)
#     --bin <sizes>            Bin sizes, comma-separated (default: 50)
#     --insertR2 <size>        Insert R2 size (default: 150)
#     --cell_num <num>         Expected cell number (default: 50000)
#     --pixelSize <size>       Pixel size in microns (default: 0.5)
# =========================

# Initialize variables with empty values
sample_name=""
chip_number=""
casno=""
chemistry=""
Species=""
method=""
mode=""

# Default optional parameters (will be set later if not provided)
OPT_THREAD=""
OPT_BIN=""
OPT_INSERTR2=""
OPT_CELL_NUM=""
OPT_PIXELSIZE=""

# Parse arguments
positional_args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            cat << 'EOF'
===============================================
Celatlas Spatial Transcriptomics Pipeline
===============================================

USAGE:
  Positional mode (legacy, backward compatible):
    ./Celatlas.sh <sample_name> <chip_number> <casno> <chemistry> <Species> <method> <mode> [options]
    ./Celatlas.sh <chip_number> <casno> <chemistry> <Species> <method> <mode> [options]

  Named arguments mode (recommended for clarity):
    ./Celatlas.sh --sample <name> --chip <chip> --casno <casno> --chemistry <chem> \
                  --species <species> --method <method> --mode <mode> [options]

  Hybrid mode (recommended for flexibility):
    ./Celatlas.sh <positional args...> --thread 32 --bin 10,50,100

REQUIRED ARGUMENTS:
  --sample, -s <name>           Sample name (optional in positional mode)
  --chip, -c <chip>             Chip/slide number (e.g., ST110250_A1)
  --casno <case>                Case/project number (e.g., HE_test)
  --chemistry <version>         Chemistry version (e.g., BBV2.4)
  --species <species>           Species name
                                  • Mus_musculus (mouse)
                                  • Homo_sapiens (human)
  --method, -m <method>         Analysis method:
                                  • image  - Image-based registration
                                  • ssDNA  - Single-strand DNA probes (alias for image)
                                  • gene_expr - Gene expression only (no image)
                                  • HE     - H&E staining guided analysis
  --mode <mode>                 Pipeline mode (typically: strna)

OPTIONAL PARAMETERS:
  --thread, -t <num>            Number of CPU threads
                                  Default: auto-detect (current: $(nproc 2>/dev/null || echo "unknown"))
                                  Recommendation: 80% of available cores

  --bin, -b <sizes>             Bin sizes in micrometers (comma-separated)
                                  Default: 50
                                  Examples:
                                    --bin 50              (single bin)
                                    --bin 10,20,50,100    (multiple bins)
                                  Guide:
                                    • 10µm  - Highest resolution, ~single cell
                                    • 20µm  - High resolution
                                    • 50µm  - Balanced (recommended)
                                    • 100µm - More genes, lower resolution

  --insertR2 <bp>               Insert R2 fragment size in base pairs
                                  Default: 150
                                  Adjust based on sequencing read length

  --cell_num <num>              Expected number of cells/spots
                                  Default: 50000

  --pixelSize <µm>              Pixel size in micrometers
                                  Default: 0.5

EXAMPLES:
  # Example 1: Basic usage with defaults
  ./Celatlas.sh DEMO ST110250_A1 HE_test BBV2.4 Mus_musculus HE strna

  # Example 2: High-performance mode (64 threads)
  ./Celatlas.sh DEMO ST110250_A1 HE_test BBV2.4 Mus_musculus HE strna --thread 64

  # Example 3: Multiple bin sizes for multi-resolution analysis
  ./Celatlas.sh DEMO ST110250_A1 HE_test BBV2.4 Mus_musculus HE strna \
                --thread 32 --bin 10,20,50,100

  # Example 4: Human sample with ssDNA probes
  ./Celatlas.sh DEMO ST110250_A1 ssDNA_test BBV2.4 Homo_sapiens ssDNA strna \
                --thread 48 --bin 50

  # Example 5: Gene expression only (no image)
  ./Celatlas.sh DEMO ST110250_A1 expr_test BBV2.4 Mus_musculus gene_expr strna

  # Example 6: Full custom configuration
  ./Celatlas.sh DEMO ST110250_A1 custom_test BBV2.4 Homo_sapiens HE strna \
                --thread 48 --bin 10,20,50,100,200 --insertR2 120 \
                --cell_num 80000 --pixelSize 0.715

  # Example 7: Pure named arguments (good for automation scripts)
  ./Celatlas.sh --sample DEMO --chip ST110250_A1 --casno batch_001 \
                --chemistry BBV2.4 --species Mus_musculus --method HE \
                --mode strna --thread 32 --bin 50,100

OUTPUT LOCATION:
  Results will be saved to:
    /mnt/strna/celatlas_spatial/results/<casno>/<chip_number>/

PIPELINE STEPS:
  01. Barcode extraction
  02. Adapter trimming (cutadapt)
  03. Alignment (STAR)
  04. Feature counting
  05. UMI counting
  06. Spatial binning (generates binX folders for each bin size)
  07. Analysis & report generation

PERFORMANCE TIPS:
  • Use --thread to set optimal CPU usage (recommended: 80% of cores)
  • Multiple bins can be generated in one run: --bin 10,50,100
  • Ensure at least 100GB free disk space
  • Recommended RAM: 128GB+ for large datasets

TROUBLESHOOTING:
  • Missing files: Check file validation messages at startup
  • Disk space: Pipeline checks and warns if < 100GB available
  • Thread count: Use `nproc` to see available CPU cores
  • Method validation: Ensure method is one of: image, ssDNA, gene_expr, HE

For detailed documentation, see:
  /mnt/strna/software/celatlas_spatial-1.6.0/USAGE.md

EOF
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
        --thread|-t)
            OPT_THREAD="$2"
            shift 2
            ;;
        --bin|-b)
            OPT_BIN="$2"
            shift 2
            ;;
        --insertR2)
            OPT_INSERTR2="$2"
            shift 2
            ;;
        --cell_num)
            OPT_CELL_NUM="$2"
            shift 2
            ;;
        --pixelSize)
            OPT_PIXELSIZE="$2"
            shift 2
            ;;
        -*)
            echo "Error: Unknown option '$1'"
            echo "Run '$0 --help' for usage information"
            exit 1
            ;;
        *)
            positional_args+=("$1")
            shift
            ;;
    esac
done

# Handle positional arguments if provided
if [ ${#positional_args[@]} -eq 7 ]; then
    sample_name="${positional_args[0]}"
    chip_number="${positional_args[1]}"
    casno="${positional_args[2]}"
    chemistry="${positional_args[3]}"
    Species="${positional_args[4]}"
    method="${positional_args[5]}"
    mode="${positional_args[6]}"
elif [ ${#positional_args[@]} -eq 6 ]; then
    sample_name=""
    chip_number="${positional_args[0]}"
    casno="${positional_args[1]}"
    chemistry="${positional_args[2]}"
    Species="${positional_args[3]}"
    method="${positional_args[4]}"
    mode="${positional_args[5]}"
elif [ ${#positional_args[@]} -gt 0 ]; then
    echo "Error: Invalid number of positional arguments (${#positional_args[@]})"
    echo "Expected 6 or 7 positional arguments, or use named arguments"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

# Validate required parameters
if [[ -z "$chip_number" || -z "$casno" || -z "$chemistry" || -z "$Species" || -z "$method" || -z "$mode" ]]; then
    echo "Error: Missing required parameters"
    echo "Run '$0 --help' for usage information"
    exit 1
fi

sample="$chip_number"

export CELATLAS_SAMPLE_NAME="$sample_name"
export CELATLAS_CHIP_NUMBER="$chip_number"

# Validate method
if [[ "$method" != "image" && "$method" != "ssDNA" && "$method" != "gene_expr" && "$method" != "HE" ]]; then
    echo "Error: Invalid method '$method'. Must be 'image', 'ssDNA', 'gene_expr', or 'HE'."
    exit 1
fi

# Map ssDNA to image for backward compatibility
if [[ "$method" == "ssDNA" ]]; then
    echo "Note: Using ssDNA mode (equivalent to image mode)"
    method="image"
fi

# =========================
# Threads & math libs (avoid numexpr cap)
# =========================
if command -v nproc >/dev/null 2>&1; then
  CPU_CORES=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
  CPU_CORES=$(sysctl -n hw.ncpu 2>/dev/null || echo 16)
else
  CPU_CORES=16
fi

# Use user-provided thread count, or fall back to environment variable, or auto-detect
if [[ -n "$OPT_THREAD" ]]; then
  thread="$OPT_THREAD"
  echo "Using user-specified thread count: $thread"
else
  thread=${THREADS:-$CPU_CORES}
  echo "Using auto-detected thread count: $thread (CPU cores: $CPU_CORES)"
fi

# Limit BLAS library threads to avoid segmentation faults
# OpenBLAS has a compiled limit (~64) but t-SNE segfaults even at 32
# Using 16 maximum for BLAS operations to ensure stability
if [ "$thread" -gt 16 ]; then
  blas_threads=16
  echo "Note: Limiting BLAS threads to 16 (user thread count: $thread)"
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

echo "[Thread Configuration]"
echo "  CPU cores detected: $CPU_CORES"
echo "  Total threads: $thread"
echo "  BLAS threads: $blas_threads (OPENBLAS_NUM_THREADS)"

# # =========================
# # Environment
# # =========================
# if [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
#   source "$HOME/anaconda3/etc/profile.d/conda.sh"
# fi
# conda activate spatial_clean || { echo "ERROR: activate env 'spatial_clean' failed"; exit 1; }

# =========================
# Config (with user overrides)
# =========================
bin=${OPT_BIN:-50}
pixelSize=${OPT_PIXELSIZE:-0.5}
insertR2=${OPT_INSERTR2:-150}
cell_num=${OPT_CELL_NUM:-50000}
feature_type="gene"
chemistryPattern="C4L15C4L15C4U10T18"
MAX_PARALLEL_FILES=${MAX_PARALLEL_FILES:-3}

# Display configuration
echo "==============================================="
echo "Celatlas Spatial Pipeline Configuration"
echo "==============================================="
echo "Sample Name:      $sample_name"
echo "Chip Number:      $chip_number"
echo "Case Number:      $casno"
echo "Chemistry:        $chemistry"
echo "Species:          $Species"
echo "Method:           $method"
echo "Mode:             $mode"
echo "-----------------------------------------------"
echo "Thread Count:     $thread"
echo "Bin Sizes:        $bin"
echo "Insert R2:        $insertR2 bp"
echo "Pixel Size:       $pixelSize µm"
echo "Expected Cells:   $cell_num"
echo "==============================================="

# Paths
workspace_dir=${CELATLAS_WORKSPACE:-/mnt/strna/celatlas_spatial}
segImageDir="${workspace_dir}/binSegment/${sample}"
image_dir="${workspace_dir}/images"
fastq_dir="${workspace_dir}/fastq/$chemistry"
mask_dir="${workspace_dir}/ST_mask"
reference_dir="${workspace_dir}/reference"
src_dir="${workspace_dir}/src"
rawdata_dir="${workspace_dir}/rawdata/$sample"
sampledir="${workspace_dir}/results/${casno}/${sample}"
genomeDir="${reference_dir}/${Species}"

# Prepare & Log
mkdir -p "${rawdata_dir}" "${sampledir}/01.barcode" "${sampledir}"
log_file="${sampledir}/pipeline.log"
exec > >(tee -a "$log_file") 2>&1
echo "Starting pipeline at $(date)"
echo "sample=${sample} casno=${casno} chemistry=${chemistry} Species=${Species} method=${method} mode=${mode} threads=${thread}"

# Validation and file preparation
echo ""
echo "[Validation] Checking required files for method=${method}, mode=${mode}..."

if [[ "$mode" == "scrna" ]]; then
    required_files=()
elif [[ "$method" == "image" ]]; then
    required_files=(
        "${image_dir}/${sample}.tif"
        "${mask_dir}/${sample}.barcodeToPos.h5"
        "${mask_dir}/${sample}_FilterBarcodes.csv"
        "${mask_dir}/${sample}_tissue_bbox.csv"
    )
    if [[ ! -f "${image_dir}/${sample}.tif" ]]; then
        echo ""
        echo "ERROR: image/ssDNA mode requires tissue image file!"
        echo "       Missing: ${image_dir}/${sample}.tif"
        echo ""
        exit 1
    fi
elif [[ "$method" == "gene_expr" ]]; then
    required_files=(
        "${mask_dir}/${sample}.barcodeToPos.h5"
        "${mask_dir}/${sample}_FilterBarcodes.csv"
        "${mask_dir}/${sample}_tissue_bbox.csv"
    )
elif [[ "$method" == "HE" ]]; then
    required_files=(
        "${mask_dir}/${sample}.barcodeToPos.h5"
        "${mask_dir}/${sample}_FilterBarcodes.csv"
        "${mask_dir}/${sample}_tissue_bbox.csv"
    )
    he_found=false
    for ext in tif png jpg jpeg; do
        if [[ -f "${image_dir}/${sample}_he.${ext}" ]]; then
            he_found=true
            break
        fi
    done
    if [[ "$he_found" == false ]]; then
        echo ""
        echo "ERROR: HE mode requires H&E staining image!"
        echo "       Missing: ${image_dir}/${sample}_he.(tif|png|jpg|jpeg)"
        echo "       Hint: Use 'gene_expr' mode if you don't have H&E image."
        echo ""
        exit 1
    fi
fi

if [[ "$mode" != "scrna" ]]; then
  for f in "${required_files[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo ""
        echo "ERROR: Required file not found: $f"
        echo ""
        exit 1
    fi
    cp -n "$f" "${rawdata_dir}/"
  done

  if [[ "$method" != "HE" ]]; then
    he_copied=false
    for ext in tif png jpg jpeg; do
      if [[ -f "${image_dir}/${sample}_he.${ext}" ]]; then
        echo "[Setup] Found optional HE image: ${sample}_he.${ext}"
        cp -n "${image_dir}/${sample}_he.${ext}" "${rawdata_dir}/"
        he_copied=true
        break
      fi
    done
  else
    for ext in tif png jpg jpeg; do
      if [[ -f "${image_dir}/${sample}_he.${ext}" ]]; then
        echo "[Setup] Copying HE image: ${sample}_he.${ext}"
        cp -n "${image_dir}/${sample}_he.${ext}" "${rawdata_dir}/"
        break
      fi
    done
  fi
fi

echo "[Validation] All required files validated successfully."
echo ""

# Reference check
[[ -d "${genomeDir}" ]] || { echo "ERROR: missing genomeDir ${genomeDir}"; exit 1; }

ulimit -n 10240

# =========================
# Detect FASTQ files (multi-lane, fold format, or simple format)
# =========================
fq1_files=""
fq2_files=""
files_found=false
sample_fq1=""

echo "Detecting FASTQ files for sample ${sample}..."

# Priority 1: Check for multi-lane sequencing format (recommended)
# Pattern: ${sample}_S*_L*_R1_*.fastq.gz
shopt -s nullglob
tenx_r1_files=( "${fastq_dir}/${sample}"_S*_L*_R1_*.fastq.gz "${fastq_dir}/${sample}"_S*_L*_R1_*.fq.gz )
tenx_r2_files=( "${fastq_dir}/${sample}"_S*_L*_R2_*.fastq.gz "${fastq_dir}/${sample}"_S*_L*_R2_*.fq.gz )
shopt -u nullglob

if [[ ${#tenx_r1_files[@]} -gt 0 && ${#tenx_r2_files[@]} -gt 0 ]]; then
  echo "✓ Detected multi-lane sequencing data (${#tenx_r1_files[@]} lane(s))"

  # Sort files by lane number
  IFS=$'\n' tenx_r1_sorted=($(sort <<<"${tenx_r1_files[*]}"))
  IFS=$'\n' tenx_r2_sorted=($(sort <<<"${tenx_r2_files[*]}"))
  unset IFS

  # Build comma-separated list
  fq1_files=$(IFS=,; echo "${tenx_r1_sorted[*]}")
  fq2_files=$(IFS=,; echo "${tenx_r2_sorted[*]}")
  sample_fq1="${tenx_r1_sorted[0]}"
  files_found=true

  echo "  R1 files: ${fq1_files}"
  echo "  R2 files: ${fq2_files}"
fi

# Priority 2: Check for fold files (legacy format)
# Pattern: ${sample}_fold1_1.fq.gz, ${sample}_fold2_1.fq.gz, ...
if [[ "$files_found" == false ]]; then
  for fold in fold1 fold2 fold3 fold4 fold5; do
    if [[ -f "${fastq_dir}/${sample}_${fold}_1.fq.gz" && -f "${fastq_dir}/${sample}_${fold}_2.fq.gz" ]]; then
      fq1_files="${fq1_files:+${fq1_files},}${fastq_dir}/${sample}_${fold}_1.fq.gz"
      fq2_files="${fq2_files:+${fq2_files},}${fastq_dir}/${sample}_${fold}_2.fq.gz"
      files_found=true
      # Use first fold file for sample step
      if [[ -z "$sample_fq1" ]]; then
        sample_fq1="${fastq_dir}/${sample}_${fold}_1.fq.gz"
      fi
    fi
  done

  if [[ "$files_found" == true ]]; then
    echo "✓ Detected multi-fold sequencing data"
    echo "  R1 files: ${fq1_files}"
    echo "  R2 files: ${fq2_files}"
  fi
fi

# Priority 3: Check for single files (simple format)
# Pattern: ${sample}_1.fq.gz, ${sample}_2.fq.gz
if [[ "$files_found" == false ]]; then
  fq1_files="${fastq_dir}/${sample}_1.fq.gz"
  fq2_files="${fastq_dir}/${sample}_2.fq.gz"
  sample_fq1="${fastq_dir}/${sample}_1.fq.gz"

  # Validate single files exist
  if [[ -f "$fq1_files" && -f "$fq2_files" ]]; then
    echo "✓ Detected single sequencing data (1 file pair)"
    echo "  R1 file: ${fq1_files}"
    echo "  R2 file: ${fq2_files}"
    files_found=true
  else
    echo "ERROR: No FASTQ files found for sample ${sample}"
    echo "  Tried multi-lane format: ${fastq_dir}/${sample}_S*_L*_R1_*.fastq.gz"
    echo "  Tried multi-fold format: ${fastq_dir}/${sample}_fold*_1.fq.gz"
    echo "  Tried single file format: ${fastq_dir}/${sample}_1.fq.gz"
    exit 1
  fi
fi

# =========================
# Validate FASTQ file integrity
# =========================
echo ""
echo "[Validation] Verifying FASTQ file integrity..."

# Parse all R1 and R2 files from comma-separated list
IFS=',' read -ra R1_ARRAY <<< "$fq1_files"
IFS=',' read -ra R2_ARRAY <<< "$fq2_files"

# Check each R1 file
for fq in "${R1_ARRAY[@]}"; do
    if [[ ! -r "$fq" ]]; then
        echo ""
        echo "ERROR: Cannot read FASTQ file: $fq"
        echo "       Please check file permissions."
        echo ""
        exit 1
    fi
    # Check if file is empty or corrupted
    if [[ ! -s "$fq" ]]; then
        echo ""
        echo "ERROR: FASTQ file is empty: $fq"
        echo ""
        exit 1
    fi
    echo "  ✓ R1: $(basename $fq) ($(du -h $fq | cut -f1))"
done

# Check each R2 file
for fq in "${R2_ARRAY[@]}"; do
    if [[ ! -r "$fq" ]]; then
        echo ""
        echo "ERROR: Cannot read FASTQ file: $fq"
        echo "       Please check file permissions."
        echo ""
        exit 1
    fi
    if [[ ! -s "$fq" ]]; then
        echo ""
        echo "ERROR: FASTQ file is empty: $fq"
        echo ""
        exit 1
    fi
    echo "  ✓ R2: $(basename $fq) ($(du -h $fq | cut -f1))"
done

echo "[Validation] All FASTQ files are readable and non-empty."
echo ""

# =========================
# Check disk space
# =========================
echo "[Disk Space] Checking available storage..."

# Get the output directory's filesystem
output_fs=$(df -h "${sampledir}" 2>/dev/null | awk 'NR==2 {print $1}')
available_space_gb=$(df -BG "${sampledir}" 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//')
available_space_human=$(df -h "${sampledir}" 2>/dev/null | awk 'NR==2 {print $4}')

# Minimum required space (GB)
required_space=100

if [[ -z "$available_space_gb" ]]; then
    echo "WARNING: Cannot determine available disk space."
    echo "         Please ensure sufficient storage is available."
else
    echo "  Filesystem: ${output_fs}"
    echo "  Available:  ${available_space_human}"
    echo "  Required:   At least ${required_space}GB recommended"

    if [[ $available_space_gb -lt $required_space ]]; then
        echo ""
        echo "WARNING: Low disk space detected!"
        echo "         Available: ${available_space_human}"
        echo "         Recommended: at least ${required_space}GB"
        echo "         Pipeline may fail if disk fills up."
        echo ""
        read -t 10 -p "Continue anyway? (y/N): " confirm || confirm="N"
        if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
            echo "Pipeline aborted by user."
            exit 1
        fi
    else
        echo "  ✓ Sufficient disk space available."
    fi
fi
echo ""

# =========================
# 00.sample
# =========================
run_command celatlas_spatial rna sample --outdir "${sampledir}/00.sample" --sample "${sample}" \
  --thread "${thread}" --chemistry "${chemistry}" --fq1 "${sample_fq1}"

# =========================
# 01.barcode
# =========================
# Determine barcode --mode for CLI (scrna|strna)
if [[ "$mode" == "scrna" ]]; then
  barcode_mode="scrna"
else
  barcode_mode="strna"
fi

# Whitelist uses the staged file in rawdata_dir for spatial
whitelist_param=""
if [[ "$barcode_mode" == "strna" ]]; then
  whitelist_param="--whitelist ${rawdata_dir}/${sample}.barcodeToPos.h5"
fi

# Determine if we need --max_parallel_files (for multiple files)
# Count the number of comma-separated files
num_files=$(echo "${fq1_files}" | tr ',' '\n' | wc -l)

if [[ $num_files -gt 1 ]]; then
  # Multiple files: use --max_parallel_files
  run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
    ${whitelist_param} --mode "${barcode_mode}" --lowNum 2 --gzip --output_R1 --resume \
    --max_parallel_files "${MAX_PARALLEL_FILES}" \
    --fq1 "${fq1_files}" --fq2 "${fq2_files}"
else
  # Single file: no need for --max_parallel_files
  run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
    ${whitelist_param} --mode "${barcode_mode}" --lowNum 2 --gzip --output_R1 --resume \
    --fq1 "${fq1_files}" --fq2 "${fq2_files}"
fi

barcode_r1="${sampledir}/01.barcode/${sample}_1.fq.gz"
barcode_r2="${sampledir}/01.barcode/${sample}_2.fq.gz"
[[ -f "$barcode_r1" && -f "$barcode_r2" ]] || { echo "ERROR: barcode outputs missing"; exit 1; }

# =========================
# 02.cutadapt  (relaxed, as in celatlas_new.sh)
# =========================
run_command celatlas_spatial rna cutadapt --outdir "${sampledir}/02.cutadapt" --sample "${sample}" \
  --thread "${thread}" --minimum_length 20 --nextseq_trim 20 --gzip --overlap 10 --insert "${insertR2}" \
  --fq "${barcode_r2}"

cutadapt_r2="${sampledir}/02.cutadapt/${sample}_clean_2.fq.gz"
[[ -f "$cutadapt_r2" ]] || { echo "ERROR: cutadapt output missing: $cutadapt_r2"; exit 1; }

# CLEAN: remove 01.barcode fastqs AFTER cutadapt succeeded
if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
  safe_rm "$barcode_r1"
  safe_rm "$barcode_r2"
fi

# =========================
# 03.STAR (adaptive parameters based on mode)
# =========================
# For scRNA mode: use more relaxed parameters to improve mapping rate
# For spatial mode: use stricter parameters for better specificity
if [[ "$mode" == "scrna" ]]; then
  # scRNA-seq: relaxed parameters to accommodate shorter/lower-quality reads
  star_match_min=30
  star_match_ratio=0.2
  star_score_ratio=0.2
  echo "[STAR] Using relaxed parameters for scRNA mode (min_match=${star_match_min}bp)"
else
  # Spatial: stricter parameters for better quality
  star_match_min=30
  star_match_ratio=0.2
  star_score_ratio=0.2
  echo "[STAR] Using strict parameters for spatial mode (min_match=${star_match_min}bp)"
fi

run_command celatlas_spatial rna star --outdir "${sampledir}/03.star" --sample "${sample}" \
  --thread "${thread}" --genomeDir "${genomeDir}" \
  --outFilterMatchNmin ${star_match_min} --outFilterMultimapNmax 1 \
  --STAR_param "--outFilterMatchNminOverLread ${star_match_ratio} --outFilterScoreMinOverLread ${star_score_ratio}" \
  --starMem 30 \
  --fq "${cutadapt_r2}"

# CLEAN: remove 02.cutadapt cleaned R2 AFTER STAR succeeded
if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
  safe_rm "$cutadapt_r2"
fi

# =========================
# 04.featureCounts
# =========================
run_command celatlas_spatial rna featureCounts --outdir "${sampledir}/04.featureCounts" --sample "${sample}" \
  --thread "${thread}" --gtf_type "${feature_type}" --genomeDir "${genomeDir}" --featureCounts_param '-s 1' \
  --input "${sampledir}/03.star/${sample}_Aligned.sortedByCoord.out.bam"

# CLEAN: remove STAR BAMs AFTER featureCounts succeeded
if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
  safe_glob_rm "${sampledir}/03.star/*.bam"
fi

# =========================
# 05.count
# =========================
run_command celatlas_spatial rna count --outdir "${sampledir}/05.count" --sample "${sample}" \
  --thread "${thread}" --genomeDir "${genomeDir}" --expected_cell_num "${cell_num}" \
  --cell_calling_method auto --bam "${sampledir}/04.featureCounts/${sample}_nameSorted.bam" --force_cell_num None

# CLEAN: remove featureCounts BAMs AFTER count succeeded
if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
  safe_glob_rm "${sampledir}/04.featureCounts/*.bam"
fi

# =========================
# 06.binSegment (branch)
# =========================
run_image_branch() {
  # Image segmentation mode: requires tissue image (auto-detected or specified)
  # Python code will auto-detect HE images with pattern: {sample}_he.(png|jpg|jpeg)
  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.binSegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${rawdata_dir}" \
    --segment --tif "${rawdata_dir}/${sample}.tif" --bs_out "${segImageDir}" --model "${src_dir}/swin_tiny.pth" --method "image" \
    --count --count_detail "${sampledir}/05.count/${sample}_count_detail.txt"
}

run_gene_expr_branch() {
  # Gene expression segmentation mode: uses UMI counts for tissue detection
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'

  echo "[Gene Expression Segmentation Mode]"
  echo "  gem_bin_size: ${gem_bin_size} microns (aggregation window)"
  echo "  umi_min_threshold: ${umi_threshold} (bottom noise pre-filtering)"
  echo "  enhance_params: ${enhance_params}"
  echo "  HE image: NO (pure gene_expr mode)"

  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.binSegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${rawdata_dir}" \
    --segment --model "${src_dir}/swin_tiny.pth" --method "gene_expr" --count \
    --count_detail "${sampledir}/05.count/${sample}_count_detail.txt" \
    --gem-bin-size "${gem_bin_size}" \
    --umi-min-threshold "${umi_threshold}" \
    --enhance-params "${enhance_params}"
}

run_he_branch() {
  # HE mode: Gene expression + HE image registration
  # Requires HE image with pattern: {sample}_he.(tif|png|jpg|jpeg)
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'
  local registration_type="${REGISTRATION_TYPE:-affine}"

  # Find HE image
  he_image_found=""
  for ext in tif png jpg jpeg; do
    if [[ -f "${rawdata_dir}/${sample}_he.${ext}" ]]; then
      he_image_found="${rawdata_dir}/${sample}_he.${ext}"
      break
    fi
  done

  if [[ -z "$he_image_found" ]]; then
    echo "[HE Mode] ERROR: No HE image found!"
    echo "          Expected: ${rawdata_dir}/${sample}_he.(tif|png|jpg|jpeg)"
    exit 1
  fi

  echo "[HE Mode: Gene Expression + HE Registration]"
  echo "  gem_bin_size: ${gem_bin_size} microns (aggregation window)"
  echo "  umi_min_threshold: ${umi_threshold} (bottom noise pre-filtering)"
  echo "  enhance_params: ${enhance_params}"
  echo "  registration_type: ${registration_type} (HE-to-GEM alignment)"
  echo "  HE image: ${he_image_found}"

  run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.binSegment" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --pixel-size "${pixelSize}" --input "${rawdata_dir}" \
    --segment --model "${src_dir}/swin_tiny.pth" --method "HE" --count \
    --count_detail "${sampledir}/05.count/${sample}_count_detail.txt" \
    --tif "${he_image_found}" \
    --gem-bin-size "${gem_bin_size}" \
    --umi-min-threshold "${umi_threshold}" \
    --enhance-params "${enhance_params}" \
    --registration-type "${registration_type}"
}

if [[ "$mode" == "scrna" ]]; then
  echo "scRNA mode: Skipping spatial steps (binSegment)."

  # =========================
  # 06_analysis_wrapper (scRNA basic analysis)
  # =========================
  echo "Running scRNA basic analysis (QC + UMAP)..."

  matrix_dir="${sampledir}/05.count/${sample}_filtered_feature_bc_matrix"

  if [[ -d "$matrix_dir" ]]; then
    run_command celatlas_spatial rna analysis \
      --outdir "${sampledir}/06_analysis_wrapper" \
      --sample "${sample}" \
      --thread "${thread}" \
      --genomeDir "${genomeDir}" \
      --matrix_file "${matrix_dir}" \
      --assay scrna

    echo "scRNA basic analysis completed: ${sampledir}/06_analysis_wrapper"
  else
    echo "WARNING: Filtered matrix not found at ${matrix_dir}, skipping analysis."
  fi

  # =========================
  # scRNA Report Generation
  # =========================
  out_html="${sample}_scrna_analysis_report.html"
  echo "Generating scRNA analysis report..."

  # Try CLI command first (if installed)
  if command -v celatlas_scrna_report >/dev/null 2>&1; then
    run_command celatlas_scrna_report \
      "${sampledir}" "${sample}" \
      --chemistry "${chemistry}" \
      --species "${Species}" \
      --output-filename "${out_html}" \
      --verbose
    echo "scRNA analysis report generated successfully: ${sampledir}/${out_html}"
  else
    # Fallback to direct Python script
    # Try to find script in multiple locations
    script_locations=(
      "/home/service/anaconda3/envs/celatlas_spatial/lib/python3.9/site-packages/celatlas_spatial/tools/scrna_report_generator.py"
      "$(dirname "$0")/celatlas_spatial/tools/scrna_report_generator.py"
      "/mnt/strna/software/celatlas_spatial-1.6.0/celatlas_spatial/tools/scrna_report_generator.py"
    )

    scrna_report_generator=""
    for loc in "${script_locations[@]}"; do
      if [[ -f "$loc" ]]; then
        scrna_report_generator="$loc"
        break
      fi
    done

    if [[ -n "$scrna_report_generator" ]]; then
      run_command python3 "$scrna_report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --output-filename "${out_html}" \
        --verbose
      echo "scRNA analysis report generated successfully (python fallback): ${sampledir}/${out_html}"
    else
      echo "WARNING: scRNA report generator not found in any of the following locations:"
      for loc in "${script_locations[@]}"; do
        echo "  - $loc"
      done
      echo "Skipping report generation."
    fi
  fi

else
  if [[ "$method" == "image" ]]; then
    run_image_branch
  elif [[ "$method" == "HE" ]]; then
    run_he_branch
  else
    run_gene_expr_branch
  fi

  # =========================
  # 07.analysis
  # =========================
  echo ""
  echo "======================================"
  echo "Step 07: Spatial Analysis"
  echo "======================================"

  # Determine dimensional reduction strategy based on bin size
  # For bin10 and bin20 (high-resolution data), skip t-SNE to avoid segfault
  # t-SNE is O(n^2) complexity and prone to memory issues with large datasets
  skip_tsne_flag=""
  if [[ "$bin" == "10" || "$bin" == "20" ]]; then
      skip_tsne_flag="--skip-tsne"
      echo "[Dimensional Reduction] Bin size: ${bin}µm - Using UMAP only (skipping t-SNE)"
      echo "  Reason: High-resolution data (bin10/bin20) has many data points"
      echo "  t-SNE is prone to segmentation fault with large datasets"
      echo "  UMAP provides similar visualization with better performance"
      echo "  Note: Report generation only requires UMAP, not t-SNE"
  else
      echo "[Dimensional Reduction] Bin size: ${bin}µm - Using both t-SNE and UMAP"
      echo "  Note: If you encounter segmentation fault, use --bin 50 or higher"
  fi
  echo ""

  run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.analysis" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --square_bin_dir "${sampledir}/06.binSegment/square_bin" \
    --pixel-size "${pixelSize}" --bin "${bin}" ${skip_tsne_flag}

 # =========================
 # 08.report (CLI first, Python fallback)
 # =========================
  out_html="${sample}_spatial_analysis_report.html"
  if command -v celatlas_spatial_report >/dev/null 2>&1; then
    run_command celatlas_spatial_report \
      "${sampledir}" "${sample}" \
      --chemistry "${chemistry}" \
      --species "${Species}" \
      --output-filename "${out_html}" \
      --verbose
    echo "Spatial analysis report generated successfully: ${sampledir}/${out_html}"
  else
    report_generator="/home/service/anaconda3/envs/celatlas_spatial/lib/python3.9/site-packages/celatlas_spatial/tools/spatial_report_generator.py"
    if [[ -f "$report_generator" ]]; then
      run_command python3 "$report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --output-filename "${out_html}" \
        --verbose
      echo "Spatial analysis report generated successfully (python fallback): ${sampledir}/${out_html}"
    else
      echo "WARNING: Report generator not found. Skipping report."
    fi
  fi
fi


echo "Pipeline completed successfully at $(date)"
