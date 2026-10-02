#!/bin/bash
set -eE

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
  CELATLAS_LAST_COMMAND="$(printf '%q ' "$@")"
  echo "Running:" "$@"
  "$@"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/celatlas_env.sh"
celatlas_env_setup "$@"
source "${SCRIPT_DIR}/scripts/celatlas_runtime_paths.sh"
source "${SCRIPT_DIR}/scripts/celatlas_post_binsegment.sh"
source "${SCRIPT_DIR}/scripts/celatlas_failure_report.sh"
celatlas_post_binsegment_init_defaults
CELATLAS_WORKFLOW="FF"

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
reference_dir=""  # New: custom reference directory
mask_dir=""       # New: custom mask directory
image_dir=""      # New: custom image directory
fastq_dir=""      # New: custom FASTQ directory
fastq_name=""     # New: custom FASTQ filename prefix

# Default optional parameters (will be set later if not provided)
OPT_THREAD=""
OPT_BIN=""
OPT_INSERTR2=""
OPT_CELL_NUM=""
OPT_PIXELSIZE=""
OPT_RESUME_EXISTING=0
OPT_BBV4_STRNA_BARCODE_MISMATCH=2
OPT_GENE_MASK_FILTER=0
OPT_STAR_MATCH_MIN="${CELATLAS_STAR_MATCH_MIN:-}"
OPT_STAR_MATCH_RATIO="${CELATLAS_STAR_MATCH_RATIO:-}"
OPT_STAR_SCORE_RATIO="${CELATLAS_STAR_SCORE_RATIO:-}"
OPT_STAR_MULTIMAP="${CELATLAS_STAR_MULTIMAP:-}"
OPT_RESOLVE_MULTIGENE_UMI="${CELATLAS_RESOLVE_MULTIGENE_UMI:-0}"

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
                                  BBV4_L9 = C20L9U10T18 compatibility library
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
                                    • 10μm  - Highest resolution, ~single cell
                                    • 20μm  - High resolution
                                    • 50μm  - Balanced (recommended)
                                    • 100μm - More genes, lower resolution

  --insertR2 <bp>               Insert R2 fragment size in base pairs
                                  Default: 150
                                  Adjust based on sequencing read length

  --cell_num <num>              Expected number of cells/spots
                                  Default: 50000

  --pixelSize <μm>              Pixel size in micrometers
                                  Default: 0.5

  --resume-existing             Reuse completed upstream outputs only when the FASTQ
                                  input manifest matches the previous successful run.
                                  If FASTQ files are added or changed, old
                                  sequencing-dependent outputs are reset and rerun.

  --bbv4-strna-barcode-mismatch <0-3>
                                Allowed barcode mismatches for BBV4/BBV4_L9
                                  spatial whitelist matching only.
                                  Default: 2. Other chemistries still use 1.

  --star-match-min <bp>         STAR --outFilterMatchNmin override.
                                  BBV4/BBV4_L9 default: 20; scrna default: 30;
                                  other spatial default: 15.

  --star-multimap <num>         STAR --outFilterMultimapNmax override.
                                  BBV4/BBV4_L9 default: 10; other default: 1.

  --resolve-multigene-umi       For barcode+UMI groups assigned to multiple genes,
                                  keep the gene with highest read support and discard
                                  ties. Useful when STAR multi-mapping is enabled.

  --gene-mask-filter            Restrict gene_expr or ssDNA/image binning to the
                                  GEM/manual gene mask region.

POST-BINSEGMENT OPTIONAL MODULES:
  --enable-cavity-filter        Run cavity/tissue mask filtering after 06.segment/01.binsegment
                                  Default: disabled

  --cavity-bins <bins>          Cavity-filtered bin list
                                  Default: 10,20,50,100,Raw

  --cavity-apply                Replace 06.segment/01.binsegment/square_bin with filtered output
                                  Default: disabled; non-destructive output is kept under
                                  06.segment/01.binsegment/square_bin_cavity_filtered
  --cavity-mask-preset <name>   Cavity mask preset: default or adipose
                                  Use adipose for H&E regions with light fat-cell interiors.

  --enable-stardist-cell-segment
                                Run StarDist H&E nuclear segmentation and barcode/count assignment
                                  Default: disabled
                                  Output: 06.segment/02.cellsegment/

  --stardist-python <path>      Python executable with stardist/tensorflow installed
                                  Default: STARDIST_PYTHON from scripts/celatlas_env.sh

  --stardist-labels <path>      Reuse an existing StarDist label TIFF
  --stardist-image <path>       DAPI/fluorescence image matching that label canvas
  --stardist-tissue-bbox <path> Label canvas bbox used for barcode coordinates

  --stardist-model-dir <path>   Local StarDist model directory
                                  Optional; can also be set in configs/celatlas.env

  --stardist-prob-thresh <num>  StarDist probability threshold
                                  Default: 0.30

  --stardist-scale <num>        StarDist internal scale
                                  Default: 2.0

  --stardist-skip-counts        Run StarDist inference/assignment only, skip count_detail aggregation
  --stardist-tiled-inference    Run full-resolution overlapping-tile inference for large TIFFs
  --stardist-tile-size <n>      Disk tile size (default: 4096)
  --stardist-tile-overlap <n>   Tile context halo (default: 256)
  --stardist-tile-merge-overlap <num>
                                Cross-tile duplicate merge fraction (default: 0.50)
  --stardist-fluorescence-channel <name>
                                gray, red, green, or blue (default: gray)
  --stardist-tile-cache-dir <path>
                                Resumable tile cache directory
  --stardist-exclude-rectangle <x0,y0,x1,y1>
                                Remove scale-bar/annotation region; repeatable
  --stardist-no-label-output    Avoid writing slide-sized expanded/QC label TIFFs

PATH CUSTOMIZATION (Advanced):
  --reference_dir <path>        Custom reference genome directory
                                  Default: $WORKSPACE/reference
                                  Useful for shared reference across projects

  --mask_dir <path>             Custom mask/barcode directory
                                  Default: $WORKSPACE/ST_mask
                                  Directory containing .barcodeToPos.h5 files

  --image_dir <path>            Custom image directory
                                  Default: $WORKSPACE/images
                                  Directory containing H&E or ssDNA images

  --fastq_dir <path>            Custom FASTQ directory
                                  Default: $WORKSPACE/fastq/$chemistry
                                  Avoids duplicate data copies

  --fastq_name <prefix>         Custom FASTQ filename prefix
                                  Default: chip_number
                                  Use when FASTQ names differ from chip number

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

  # Example 8: Custom directories for shared data environment
  ./Celatlas.sh --chip ST110250_A1 --casno shared_001 \
                --chemistry BBV2.4 --species Mus_musculus --method HE \
                --mode strna \
                --reference_dir /data/shared/reference \
                --fastq_dir /data/sequencing/run001 \
                --mask_dir /data/shared/masks \
                --image_dir /data/shared/images

  # Example 9: HE mode with post-binSegment cavityFilter + StarDist cell matrix
  ./Celatlas.sh DEMO ST110250_A1 HE_test BBV2.4 Mus_musculus HE strna \
                --thread 32 --bin 50 \
                --enable-cavity-filter \
                --enable-stardist-cell-segment

OUTPUT LOCATION:
  Results will be saved to:
    $CELATLAS_RESULTS_ROOT/<casno>/<chip_number>/

PIPELINE STEPS:
  01. Barcode extraction
  02. Adapter trimming (cutadapt)
  03. Alignment (STAR)
  04. Feature counting
  05. UMI counting
  06. Spatial binning (generates binX folders for each bin size)
  06.segment/01 Optional cavityFilter (--enable-cavity-filter)
  06.segment/02 Optional StarDist cell matrix (--enable-stardist-cell-segment)
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

CONFIGURATION:
  Copy configs/celatlas.env.example to configs/celatlas.env or celatlas.local.env
  to customize workspace, reference, environment, and model paths.

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
        --reference_dir)
            reference_dir="$2"
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
        --fastq_dir)
            fastq_dir="$2"
            shift 2
            ;;
        --fastq_name)
            fastq_name="$2"
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
        --resume-existing)
            OPT_RESUME_EXISTING=1
            shift
            ;;
        --bbv4-strna-barcode-mismatch|--bbv4-strna-mismatch)
            OPT_BBV4_STRNA_BARCODE_MISMATCH="$2"
            shift 2
            ;;
        --star-match-min)
            OPT_STAR_MATCH_MIN="$2"
            shift 2
            ;;
        --star-match-ratio)
            OPT_STAR_MATCH_RATIO="$2"
            shift 2
            ;;
        --star-score-ratio)
            OPT_STAR_SCORE_RATIO="$2"
            shift 2
            ;;
        --star-multimap)
            OPT_STAR_MULTIMAP="$2"
            shift 2
            ;;
        --resolve-multigene-umi)
            OPT_RESOLVE_MULTIGENE_UMI=1
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
if [ ${#positional_args[@]} -eq 12 ]; then
    # Full mode: all parameters including custom directories
    sample_name="${positional_args[0]}"
    chip_number="${positional_args[1]}"
    casno="${positional_args[2]}"
    chemistry="${positional_args[3]}"
    Species="${positional_args[4]}"
    method="${positional_args[5]}"
    mode="${positional_args[6]}"
    reference_dir="${positional_args[7]}"
    mask_dir="${positional_args[8]}"
    image_dir="${positional_args[9]}"
    fastq_dir="${positional_args[10]}"
    fastq_name="${positional_args[11]}"
elif [ ${#positional_args[@]} -eq 11 ]; then
    # No image_dir specified (for gene_expr mode)
    sample_name="${positional_args[0]}"
    chip_number="${positional_args[1]}"
    casno="${positional_args[2]}"
    chemistry="${positional_args[3]}"
    Species="${positional_args[4]}"
    method="${positional_args[5]}"
    mode="${positional_args[6]}"
    reference_dir="${positional_args[7]}"
    mask_dir="${positional_args[8]}"
    image_dir=""
    fastq_dir="${positional_args[9]}"
    fastq_name="${positional_args[10]}"
elif [ ${#positional_args[@]} -eq 7 ]; then
    # Standard mode with sample_name
    sample_name="${positional_args[0]}"
    chip_number="${positional_args[1]}"
    casno="${positional_args[2]}"
    chemistry="${positional_args[3]}"
    Species="${positional_args[4]}"
    method="${positional_args[5]}"
    mode="${positional_args[6]}"
    reference_dir=""
    mask_dir=""
    image_dir=""
    fastq_dir=""
    fastq_name=""
elif [ ${#positional_args[@]} -eq 6 ]; then
    # Legacy mode without sample_name
    sample_name=""
    chip_number="${positional_args[0]}"
    casno="${positional_args[1]}"
    chemistry="${positional_args[2]}"
    Species="${positional_args[3]}"
    method="${positional_args[4]}"
    mode="${positional_args[5]}"
    reference_dir=""
    mask_dir=""
    image_dir=""
    fastq_dir=""
    fastq_name=""
elif [ ${#positional_args[@]} -gt 0 ]; then
    echo "Error: Invalid number of positional arguments (${#positional_args[@]})"
    echo "Expected 6, 7, 11, or 12 positional arguments, or use named arguments"
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

fastq_chemistry="$chemistry"
if [[ "$chemistry" == "BBV4_L9" ]]; then
    # Compatibility library: the FASTQ directory remains grouped under BBV4.
    fastq_chemistry="BBV4"
fi

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
celatlas_post_binsegment_finalize_defaults "$method"

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
if [[ "$chemistry" == "BBV4" ]]; then
  chemistryPattern="C20L15U10T18"
elif [[ "$chemistry" == "BBV4_L9" ]]; then
  chemistryPattern="C20L9U10T18"
fi
MAX_PARALLEL_FILES=${MAX_PARALLEL_FILES:-3}

# Paths - set defaults if not provided by user/config
workspace_dir="${CELATLAS_WORKSPACE}"

# Apply defaults for custom directories if not specified
if [[ -z "$reference_dir" ]]; then
  reference_dir="${CELATLAS_REFERENCE_DIR}"
  echo "Note: Using default reference directory: $reference_dir"
fi

if [[ -z "$mask_dir" ]]; then
  mask_dir="${CELATLAS_MASK_DIR}"
  echo "Note: Using default mask directory: $mask_dir"
fi

if [[ -z "$image_dir" ]]; then
  image_dir="${CELATLAS_IMAGE_DIR}"
  echo "Note: Using default image directory: $image_dir"
fi

if [[ -z "$fastq_dir" ]]; then
  fastq_dir="${CELATLAS_FASTQ_ROOT}/$fastq_chemistry"
  echo "Note: Using default FASTQ directory: $fastq_dir"
fi

if [[ -z "$fastq_name" ]]; then
  fastq_name="$sample"
  echo "Note: FASTQ name not provided, using chip number: $fastq_name"
fi

src_dir="${CELATLAS_SRC_DIR}"
sampledir="${CELATLAS_RESULTS_ROOT}/${casno}/${sample}"
runtime_mask_dir=""
if [[ "$mode" != "scrna" ]]; then
  runtime_mask_dir="$(celatlas_mask_dir_for_sample "${sampledir}")"
fi
genomeDir="${reference_dir}/${Species}"

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
echo "Reference Dir:    $reference_dir"
echo "Mask Dir:         $mask_dir"
echo "Image Dir:        $image_dir"
echo "FASTQ Dir:        $fastq_dir"
echo "FASTQ Name:       $fastq_name"
echo "-----------------------------------------------"
echo "Thread Count:     $thread"
echo "Bin Sizes:        $bin"
echo "Insert R2:        $insertR2 bp"
echo "Pixel Size:       $pixelSize μm"
echo "Expected Cells:   $cell_num"
echo "Resume Existing:  $OPT_RESUME_EXISTING"
echo "BBV4 strna BC MM: $OPT_BBV4_STRNA_BARCODE_MISMATCH"
echo "Resolve multi UMI:$OPT_RESOLVE_MULTIGENE_UMI"
echo "Gene Mask Filter: $OPT_GENE_MASK_FILTER"
echo "Cavity Filter:    $OPT_ENABLE_CAVITY_FILTER"
echo "StarDist Segment: $OPT_ENABLE_STARDIST_CELL_SEGMENT"
if [[ "$OPT_ENABLE_CAVITY_FILTER" == "1" ]]; then
  echo "Cavity Bins:      $OPT_CAVITY_BINS"
fi
if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" == "1" ]]; then
  echo "StarDist Model:   $OPT_STARDIST_MODEL"
  echo "StarDist p/scale: ${OPT_STARDIST_PROB_THRESH}/${OPT_STARDIST_SCALE}"
fi
echo "==============================================="

# Prepare & Log
if [[ "$mode" == "scrna" ]]; then
  mkdir -p "${sampledir}/01.barcode" "${sampledir}"
else
  mkdir -p "${runtime_mask_dir}" "${sampledir}/01.barcode" "${sampledir}"
fi
log_file="${sampledir}/pipeline.log"
exec > >(tee -a "$log_file") 2>&1
trap 'celatlas_on_pipeline_error "$?" "$LINENO" "$BASH_COMMAND"' ERR
trap 'celatlas_on_pipeline_exit "$?" "$LINENO" "$BASH_COMMAND"' EXIT
echo "Starting pipeline at $(date)"
echo "sample=${sample} casno=${casno} chemistry=${chemistry} Species=${Species} method=${method} mode=${mode} threads=${thread} reference_dir=${reference_dir} mask_dir=${mask_dir} image_dir=${image_dir} fastq_dir=${fastq_dir} fastq_name=${fastq_name} enable_cavity_filter=${OPT_ENABLE_CAVITY_FILTER} enable_stardist_cell_segment=${OPT_ENABLE_STARDIST_CELL_SEGMENT}"
echo "runtime_mask_dir=${runtime_mask_dir:-N/A}"
echo "resume_existing=${OPT_RESUME_EXISTING}"

# Validation and file preparation
echo ""
echo "[Validation] Checking required files for method=${method}, mode=${mode}..."

if [[ "$mode" == "scrna" ]]; then
    required_files=()
elif [[ "$method" == "image" ]]; then
    required_files=(
        "${image_dir}/${sample}.tif"
    )
    if [[ ! -f "${image_dir}/${sample}.tif" ]]; then
        echo ""
        echo "ERROR: image/ssDNA mode requires tissue image file!"
        echo "       Missing: ${image_dir}/${sample}.tif"
        echo ""
        exit 1
    fi
elif [[ "$method" == "gene_expr" ]]; then
    required_files=()
elif [[ "$method" == "HE" ]]; then
    required_files=()
    he_found=false
    for ext in tif tiff png jpg jpeg; do
        if [[ -f "${image_dir}/${sample}_he.${ext}" ]]; then
            he_found=true
            break
        fi
    done
    if [[ "$he_found" == false ]]; then
        echo ""
        echo "ERROR: HE mode requires H&E staining image!"
        echo "       Missing: ${image_dir}/${sample}_he.(tif|tiff|png|jpg|jpeg)"
        echo "       Hint: Use 'gene_expr' mode if you don't have H&E image."
        echo ""
        exit 1
    fi
fi

if [[ "$mode" != "scrna" ]]; then
  celatlas_stage_spatial_inputs "$sample" "$sampledir" "$mask_dir" "$image_dir" "$method"
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

echo "Detecting FASTQ files for FASTQ name: ${fastq_name}..."

# Priority 1: Check for multi-lane sequencing format (recommended)
# Pattern: ${fastq_name}_S*_L*_R1_*.fastq.gz
shopt -s nullglob
tenx_r1_files=( "${fastq_dir}/${fastq_name}"_S*_L*_R1_*.fastq.gz "${fastq_dir}/${fastq_name}"_S*_L*_R1_*.fq.gz )
tenx_r2_files=( "${fastq_dir}/${fastq_name}"_S*_L*_R2_*.fastq.gz "${fastq_dir}/${fastq_name}"_S*_L*_R2_*.fq.gz )
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
# Pattern: ${fastq_name}_fold1_1.fq.gz, ${fastq_name}_fold2_1.fq.gz, ...
if [[ "$files_found" == false ]]; then
  for fold in fold1 fold2 fold3 fold4 fold5; do
    if [[ -f "${fastq_dir}/${fastq_name}_${fold}_1.fq.gz" && -f "${fastq_dir}/${fastq_name}_${fold}_2.fq.gz" ]]; then
      fq1_files="${fq1_files:+${fq1_files},}${fastq_dir}/${fastq_name}_${fold}_1.fq.gz"
      fq2_files="${fq2_files:+${fq2_files},}${fastq_dir}/${fastq_name}_${fold}_2.fq.gz"
      files_found=true
      # Use first fold file for sample step
      if [[ -z "$sample_fq1" ]]; then
        sample_fq1="${fastq_dir}/${fastq_name}_${fold}_1.fq.gz"
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
# Pattern: ${fastq_name}_1.fq.gz, ${fastq_name}_2.fq.gz
if [[ "$files_found" == false ]]; then
  fq1_files="${fastq_dir}/${fastq_name}_1.fq.gz"
  fq2_files="${fastq_dir}/${fastq_name}_2.fq.gz"
  sample_fq1="${fastq_dir}/${fastq_name}_1.fq.gz"

  # Validate single files exist
  if [[ -f "$fq1_files" && -f "$fq2_files" ]]; then
    echo "✓ Detected single sequencing data (1 file pair)"
    echo "  R1 file: ${fq1_files}"
    echo "  R2 file: ${fq2_files}"
    files_found=true
  fi
fi

# Priority 4: Check for _R1/_R2 format (common alternative format)
# Pattern: ${fastq_name}_R1.fq.gz, ${fastq_name}_R2.fq.gz
if [[ "$files_found" == false ]]; then
  fq1_files="${fastq_dir}/${fastq_name}_R1.fq.gz"
  fq2_files="${fastq_dir}/${fastq_name}_R2.fq.gz"
  sample_fq1="${fastq_dir}/${fastq_name}_R1.fq.gz"

  # Validate single files exist
  if [[ -f "$fq1_files" && -f "$fq2_files" ]]; then
    echo "✓ Detected _R1/_R2 format sequencing data (1 file pair)"
    echo "  R1 file: ${fq1_files}"
    echo "  R2 file: ${fq2_files}"
    files_found=true
  fi
fi

# Final check: no files found
if [[ "$files_found" == false ]]; then
  echo "ERROR: No FASTQ files found for FASTQ name: ${fastq_name}"
  echo "  Tried multi-lane format: ${fastq_dir}/${fastq_name}_S*_L*_R1_*.fastq.gz"
  echo "  Tried multi-fold format: ${fastq_dir}/${fastq_name}_fold*_1.fq.gz"
  echo "  Tried simple format: ${fastq_dir}/${fastq_name}_1.fq.gz"
  echo "  Tried _R1/_R2 format: ${fastq_dir}/${fastq_name}_R1.fq.gz"
  echo ""
  echo "Please ensure FASTQ files exist in: ${fastq_dir}"
  echo "Or specify custom FASTQ directory with --fastq_dir option"
  echo "Or specify custom FASTQ name prefix with --fastq_name option"
  exit 1
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

fastq_manifest_file="${sampledir}/.fastq_inputs.tsv"
fastq_manifest_current="${sampledir}/.fastq_inputs.current.tsv"
write_fastq_manifest() {
  local output="$1"
  : > "$output"
  local fq
  for fq in "${R1_ARRAY[@]}"; do
    printf 'R1\t%s\t%s\t%s\n' "$fq" "$(stat -c %s "$fq")" "$(stat -c %Y "$fq")" >> "$output"
  done
  for fq in "${R2_ARRAY[@]}"; do
    printf 'R2\t%s\t%s\t%s\n' "$fq" "$(stat -c %s "$fq")" "$(stat -c %Y "$fq")" >> "$output"
  done
}
write_fastq_manifest "$fastq_manifest_current"

FASTQ_MANIFEST_STATUS="$(celatlas_fastq_manifest_status "$fastq_manifest_file" "$fastq_manifest_current")"
RESUME_UPSTREAM=0
ARCHIVE_FASTQ_DEPENDENT=0

case "$FASTQ_MANIFEST_STATUS" in
  same)
    if [[ "$OPT_RESUME_EXISTING" == "1" ]]; then
      RESUME_UPSTREAM=1
      echo "[Resume] FASTQ manifest unchanged; completed upstream stages may be reused."
    else
      echo "[FASTQ] FASTQ manifest unchanged; upstream stages will rerun because --resume-existing was not set."
    fi
    ;;
  added)
    echo "[FASTQ] Detected supplementary sequencing data: current manifest adds FASTQ files to the previous run."
    echo "[FASTQ] Sequencing-dependent outputs will be reset and regenerated to keep UMI deduplication correct."
    ARCHIVE_FASTQ_DEPENDENT=1
    ;;
  reordered)
    echo "[FASTQ] FASTQ manifest order changed; treating as changed for a clean rerun."
    ARCHIVE_FASTQ_DEPENDENT=1
    ;;
  changed)
    echo "[FASTQ] FASTQ manifest changed; sequencing-dependent outputs will be reset and regenerated."
    ARCHIVE_FASTQ_DEPENDENT=1
    ;;
  new)
    if celatlas_has_fastq_dependent_outputs "$sampledir" "$sample"; then
      echo "[FASTQ] Previous FASTQ manifest missing but sequencing-dependent outputs exist."
      ARCHIVE_FASTQ_DEPENDENT=1
    else
      echo "[FASTQ] No previous FASTQ manifest; running upstream stages."
    fi
    ;;
esac

if [[ "$ARCHIVE_FASTQ_DEPENDENT" == "1" ]]; then
  celatlas_reset_fastq_dependent_outputs "$sampledir" "$sample" "$FASTQ_MANIFEST_STATUS"
fi

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
if [[ "$RESUME_UPSTREAM" == "1" && -f "${sampledir}/00.sample/step.log" ]]; then
  echo "[Resume] Reusing 00.sample"
else
  run_command celatlas_spatial rna sample --outdir "${sampledir}/00.sample" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --fq1 "${sample_fq1}"
fi

# =========================
# 01.barcode
# =========================
# Determine barcode --mode for CLI (scrna|strna)
if [[ "$mode" == "scrna" ]]; then
  barcode_mode="scrna"
else
  barcode_mode="strna"
fi

# Whitelist 选择：
# - strna 模式：优先使用 .barcodeToPos.h5；缺失时使用 _FilterBarcodes.csv
#   _tissue_bbox.csv 仅作为 legacy optional hint
# - scrna + BBV4：第一阶段不传 whitelist（barcode 100% 提取）
# - scrna + 其他化学：使用 chemistry 内置 bclist（不传 --whitelist）
whitelist_param=""
if [[ "$barcode_mode" == "strna" ]]; then
  if [[ -f "${runtime_mask_dir}/${sample}.barcodeToPos.h5" ]]; then
    whitelist_param="--whitelist ${runtime_mask_dir}/${sample}.barcodeToPos.h5"
  else
    whitelist_param="--whitelist ${runtime_mask_dir}/${sample}_FilterBarcodes.csv"
    if [[ -f "${runtime_mask_dir}/${sample}_tissue_bbox.csv" ]]; then
      whitelist_param="${whitelist_param} --tissue-bbox ${runtime_mask_dir}/${sample}_tissue_bbox.csv"
    fi
  fi
elif [[ "$barcode_mode" == "scrna" && "$chemistry" =~ ^BBV4(_L9)?$ ]]; then
  echo "Note: ${chemistry} scrna 第一阶段不使用 whitelist，仅进行 linker 过滤"
  whitelist_param=""
fi

# Determine if we need --max_parallel_files (for multiple files)
# Count the number of comma-separated files
num_files=$(echo "${fq1_files}" | tr ',' '\n' | wc -l)

barcode_r1="${sampledir}/01.barcode/${sample}_1.fq.gz"
barcode_r2="${sampledir}/01.barcode/${sample}_2.fq.gz"
if [[ "$RESUME_UPSTREAM" == "1" && -f "$barcode_r1" && -f "$barcode_r2" ]]; then
  echo "[Resume] Reusing 01.barcode"
elif [[ $num_files -gt 1 ]]; then
  # Multiple files: use --max_parallel_files
  run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
    ${whitelist_param} --mode "${barcode_mode}" --lowNum 2 --gzip --output_R1 --resume \
    --bbv4-strna-barcode-mismatch "${OPT_BBV4_STRNA_BARCODE_MISMATCH}" \
    --max_parallel_files "${MAX_PARALLEL_FILES}" \
    --fq1 "${fq1_files}" --fq2 "${fq2_files}"
else
  # Single file: no need for --max_parallel_files
  run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
    ${whitelist_param} --mode "${barcode_mode}" --lowNum 2 --gzip --output_R1 --resume \
    --bbv4-strna-barcode-mismatch "${OPT_BBV4_STRNA_BARCODE_MISMATCH}" \
    --fq1 "${fq1_files}" --fq2 "${fq2_files}"
fi

[[ -f "$barcode_r1" && -f "$barcode_r2" ]] || { echo "ERROR: barcode outputs missing"; exit 1; }

# =========================
# 02.cutadapt  (relaxed, as in celatlas_new.sh)
# =========================
cutadapt_r2="${sampledir}/02.cutadapt/${sample}_clean_2.fq.gz"
if [[ "$RESUME_UPSTREAM" == "1" && -f "$cutadapt_r2" ]]; then
  echo "[Resume] Reusing 02.cutadapt"
else
  run_command celatlas_spatial rna cutadapt --outdir "${sampledir}/02.cutadapt" --sample "${sample}" \
    --thread "${thread}" --minimum_length 20 --nextseq_trim 20 --gzip --overlap 10 --insert "${insertR2}" \
    --fq "${barcode_r2}"
fi
[[ -f "$cutadapt_r2" ]] || { echo "ERROR: cutadapt output missing: $cutadapt_r2"; exit 1; }

# CLEAN: remove 01.barcode fastqs AFTER cutadapt succeeded
# if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
#   safe_rm "$barcode_r1"
#   safe_rm "$barcode_r2"
# fi

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
  star_multimap=1
  echo "[STAR] Using relaxed parameters for scRNA mode (min_match=${star_match_min}bp)"
elif [[ "$chemistry" == "BBV4" || "$chemistry" == "BBV4_L9" ]]; then
  # HD BBV4 R2 often has valid shorter genomic inserts after polyA/adapter trimming.
  star_match_min=20
  star_match_ratio=0.2
  star_score_ratio=0.2
  star_multimap=10
  echo "[STAR] Using HD BBV4 parameters (min_match=${star_match_min}bp, multimap=${star_multimap})"
else
  # Spatial: stricter parameters for better quality
  star_match_min=15
  star_match_ratio=0.2
  star_score_ratio=0.2
  star_multimap=1
  echo "[STAR] Using strict parameters for spatial mode (min_match=${star_match_min}bp)"
fi

star_match_min="${OPT_STAR_MATCH_MIN:-$star_match_min}"
star_match_ratio="${OPT_STAR_MATCH_RATIO:-$star_match_ratio}"
star_score_ratio="${OPT_STAR_SCORE_RATIO:-$star_score_ratio}"
star_multimap="${OPT_STAR_MULTIMAP:-$star_multimap}"
echo "[STAR] Effective parameters: min_match=${star_match_min}, multimap=${star_multimap}, match_ratio=${star_match_ratio}, score_ratio=${star_score_ratio}"

star_bam="${sampledir}/03.star/${sample}_Aligned.sortedByCoord.out.bam"
if [[ "$RESUME_UPSTREAM" == "1" && -f "$star_bam" ]]; then
  echo "[Resume] Reusing 03.star"
else
  run_command celatlas_spatial rna star --outdir "${sampledir}/03.star" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" \
    --outFilterMatchNmin "${star_match_min}" --outFilterMultimapNmax "${star_multimap}" \
    --STAR_param "--outFilterMatchNminOverLread ${star_match_ratio} --outFilterScoreMinOverLread ${star_score_ratio}" \
    --starMem 30 \
    --fq "${cutadapt_r2}"
fi

# # CLEAN: remove 02.cutadapt cleaned R2 AFTER STAR succeeded
# if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
#   safe_rm "$cutadapt_r2"
# fi

# =========================
# 04.featureCounts
# =========================
featurecounts_bam="${sampledir}/04.featureCounts/${sample}_nameSorted.bam"
if [[ "$RESUME_UPSTREAM" == "1" && -f "$featurecounts_bam" ]]; then
  echo "[Resume] Reusing 04.featureCounts"
else
  run_command celatlas_spatial rna featureCounts --outdir "${sampledir}/04.featureCounts" --sample "${sample}" \
    --thread "${thread}" --gtf_type "${feature_type}" --genomeDir "${genomeDir}" --featureCounts_param '-s 1' \
    --input "${star_bam}"
fi

# # CLEAN: remove STAR BAMs AFTER featureCounts succeeded
# if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
#   safe_glob_rm "${sampledir}/03.star/*.bam"
# fi

# =========================
# 05.count
# =========================
count_detail="${sampledir}/05.count/${sample}_count_detail.txt"
if [[ "$RESUME_UPSTREAM" == "1" && -f "$count_detail" ]]; then
  echo "[Resume] Reusing 05.count"
else
  count_extra_args=()
  if [[ "$OPT_RESOLVE_MULTIGENE_UMI" == "1" || "$OPT_RESOLVE_MULTIGENE_UMI" == "true" || "$OPT_RESOLVE_MULTIGENE_UMI" == "TRUE" ]]; then
    count_extra_args+=(--resolve-multigene-umi)
  fi
  run_command celatlas_spatial rna count --outdir "${sampledir}/05.count" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --expected_cell_num "${cell_num}" \
    --cell_calling_method auto --bam "${featurecounts_bam}" --force_cell_num None \
    "${count_extra_args[@]}"
  cp -f "$fastq_manifest_current" "$fastq_manifest_file"
fi

# # CLEAN: remove featureCounts BAMs AFTER count succeeded
# if [[ "$CLEAN_INTERMEDIATE" == "1" ]]; then
#   safe_glob_rm "${sampledir}/04.featureCounts/*.bam"
# fi

# =========================
# 06.segment/01.binsegment (branch)
# =========================
run_image_branch() {
  # Image segmentation mode: requires tissue image (auto-detected or specified)
  # Python code will auto-detect HE images with pattern: {sample}_he.(png|jpg|jpeg)
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
  # Gene expression segmentation mode: uses UMI counts for tissue detection
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'

  echo "[Gene Expression Segmentation Mode]"
  echo "  gem_bin_size: ${gem_bin_size} microns (aggregation window)"
  echo "  umi_min_threshold: ${umi_threshold} (bottom noise pre-filtering)"
  echo "  enhance_params: ${enhance_params}"
  echo "  HE image: NO (pure gene_expr mode)"

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
  # HE mode: Gene expression + HE image registration
  # Requires HE image with pattern: {sample}_he.(tif|png|jpg|jpeg)
  local gem_bin_size="${GEM_BIN_SIZE:-20}"
  local umi_threshold="${UMI_MIN_THRESHOLD:-30}"
  local enhance_params='{"p_low":5, "p_high":95, "suppress_noise":true, "noise_threshold":"auto"}'
  local registration_type="${REGISTRATION_TYPE:-affine}"

  # Find HE image
  he_image_found="$(celatlas_find_runtime_he_image "$sample" "$sampledir" "$image_dir" || true)"

  if [[ -z "$he_image_found" ]]; then
    echo "[HE Mode] ERROR: No HE image found!"
    echo "          Expected: ${runtime_mask_dir}/${sample}_he.(tif|tiff|png|jpg|jpeg)"
    exit 1
  fi

  echo "[HE Mode: Gene Expression + HE Registration]"
  echo "  gem_bin_size: ${gem_bin_size} microns (aggregation window)"
  echo "  umi_min_threshold: ${umi_threshold} (bottom noise pre-filtering)"
  echo "  enhance_params: ${enhance_params}"
  echo "  registration_type: ${registration_type} (HE-to-GEM alignment)"
  echo "  HE image: ${he_image_found}"

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

run_post_binsegment_optional_modules() {
  celatlas_run_post_binsegment_optional_modules
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
    # Fallback to direct Python script from the active checkout or config.
    scrna_report_generator="${CELATLAS_SCRNA_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/scrna_report_generator.py}"

    if [[ -f "$scrna_report_generator" ]]; then
      run_command python3 "$scrna_report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --output-filename "${out_html}" \
        --verbose
      echo "scRNA analysis report generated successfully (python fallback): ${sampledir}/${out_html}"
    else
      echo "WARNING: scRNA report generator not found: ${scrna_report_generator}"
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

  celatlas_sync_gem_expression_to_mask_dir "$sampledir"

  run_post_binsegment_optional_modules

  # =========================
  # 07.outs
  # =========================
  echo ""
  echo "======================================"
  echo "Step 07: Outs Analysis"
  echo "======================================"

  # Determine dimensional reduction strategy based on bin size
  # For bin10 and bin20 (high-resolution data), skip t-SNE to avoid segfault
  # t-SNE is O(n^2) complexity and prone to memory issues with large datasets
  skip_tsne_flag=""
  if [[ "$bin" == "10" || "$bin" == "20" ]]; then
      skip_tsne_flag="--skip-tsne"
      echo "[Dimensional Reduction] Bin size: ${bin}μm - Using UMAP only (skipping t-SNE)"
      echo "  Reason: High-resolution data (bin10/bin20) has many data points"
      echo "  t-SNE is prone to segmentation fault with large datasets"
      echo "  UMAP provides similar visualization with better performance"
      echo "  Note: Report generation only requires UMAP, not t-SNE"
  else
      echo "[Dimensional Reduction] Bin size: ${bin}μm - Using both t-SNE and UMAP"
      echo "  Note: If you encounter segmentation fault, use --bin 50 or higher"
  fi
  echo ""

  run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.outs/binned_outputs" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${genomeDir}" --square_bin_dir "${sampledir}/06.segment/01.binsegment/square_bin" \
    --pixel-size "${pixelSize}" --bin "${bin}" ${skip_tsne_flag}

  cell_matrix_dir="${sampledir}/06.segment/02.cellsegment/cell_matrix"
  if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" == "1" && -d "${cell_matrix_dir}" && -f "${cell_matrix_dir}/matrix.mtx.gz" ]]; then
    echo ""
    echo "======================================"
    echo "Step 07: Cell-Segmented Analysis"
    echo "======================================"
    run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.outs/cellsegmented_outputs" --sample "${sample}" \
      --thread "${thread}" --genomeDir "${genomeDir}" --assay cell \
      --cell-seg-dir "${sampledir}/06.segment/02.cellsegment" \
      --skip-tsne
  fi

 # =========================
 # 08.report (CLI first, Python fallback)
 # =========================
  out_html="${sample}_spatial_analysis_report.html"
  report_cellsegment_args=()
  if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" != "1" ]]; then
    report_cellsegment_args+=(--exclude-cell-segmentation)
  fi
  if command -v celatlas_spatial_report >/dev/null 2>&1; then
    run_command celatlas_spatial_report \
      "${sampledir}" "${sample}" \
      --chemistry "${chemistry}" \
      --species "${Species}" \
      --method "${method}" \
      --output-filename "${out_html}" \
      --verbose \
      "${report_cellsegment_args[@]}"
    echo "Spatial analysis report generated successfully: ${sampledir}/${out_html}"
  else
    report_generator="${CELATLAS_SPATIAL_REPORT_GENERATOR:-${SCRIPT_DIR}/celatlas_spatial/tools/spatial_report_generator.py}"
    if [[ -f "$report_generator" ]]; then
      run_command python3 "$report_generator" \
        "${sampledir}" \
        "${sample}" \
        --chemistry "${chemistry}" \
        --species "${Species}" \
        --method "${method}" \
        --output-filename "${out_html}" \
        --verbose \
        "${report_cellsegment_args[@]}"
      echo "Spatial analysis report generated successfully (python fallback): ${sampledir}/${out_html}"
    else
      echo "WARNING: Report generator not found: ${report_generator}. Skipping report."
    fi
  fi
fi


echo "Pipeline completed successfully at $(date)"
