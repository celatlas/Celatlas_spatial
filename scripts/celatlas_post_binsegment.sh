#!/usr/bin/env bash

# Shared post-binSegment module for the Celatlas v1.8 shell entrypoints.
# The caller must define sampledir, sample, method, and run_command.

celatlas_post_binsegment_init_defaults() {
  if [[ -z "${OPT_ENABLE_CAVITY_FILTER+x}" ]]; then
    OPT_ENABLE_CAVITY_FILTER="${ENABLE_CAVITY_FILTER:-0}"
  fi
  if [[ -z "${OPT_ENABLE_STARDIST_CELL_SEGMENT+x}" ]]; then
    OPT_ENABLE_STARDIST_CELL_SEGMENT="${ENABLE_STARDIST_CELL_SEGMENT:-${ENABLE_STARDIST:-0}}"
  fi
  if [[ -z "${OPT_CAVITY_BINS+x}" ]]; then
    OPT_CAVITY_BINS="${CAVITY_BINS:-10,20,50,100,Raw}"
  fi
  if [[ -z "${OPT_CAVITY_APPLY+x}" ]]; then
    OPT_CAVITY_APPLY="${CAVITY_APPLY:-0}"
  fi
  if [[ -z "${OPT_CAVITY_DRY_RUN+x}" ]]; then
    OPT_CAVITY_DRY_RUN="${CAVITY_DRY_RUN:-0}"
  fi
  if [[ -z "${OPT_CAVITY_SKIP_QC_IMAGES+x}" ]]; then
    OPT_CAVITY_SKIP_QC_IMAGES="${CAVITY_SKIP_QC_IMAGES:-0}"
  fi
  if [[ -z "${OPT_CAVITY_MASK_PRESET+x}" ]]; then
    OPT_CAVITY_MASK_PRESET="${CAVITY_MASK_PRESET:-${MASK_PRESET:-default}}"
  fi
  if [[ -z "${OPT_STARDIST_PYTHON+x}" ]]; then
    OPT_STARDIST_PYTHON="${STARDIST_PYTHON:-}"
  fi
  if [[ -z "${OPT_STARDIST_LABELS+x}" ]]; then
    OPT_STARDIST_LABELS="${STARDIST_LABELS:-}"
  fi
  if [[ -z "${OPT_STARDIST_IMAGE+x}" ]]; then
    OPT_STARDIST_IMAGE="${STARDIST_IMAGE:-}"
  fi
  if [[ -z "${OPT_STARDIST_TISSUE_BBOX+x}" ]]; then
    OPT_STARDIST_TISSUE_BBOX="${STARDIST_TISSUE_BBOX:-}"
  fi
  if [[ -z "${OPT_STARDIST_MODEL+x}" ]]; then
    OPT_STARDIST_MODEL="${STARDIST_MODEL:-}"
  fi
  if [[ -z "${OPT_STARDIST_MODEL_DIR+x}" ]]; then
    OPT_STARDIST_MODEL_DIR="${STARDIST_MODEL_DIR:-}"
  fi
  if [[ -z "${OPT_STARDIST_PROB_THRESH+x}" ]]; then
    OPT_STARDIST_PROB_THRESH="${STARDIST_PROB_THRESH:-0.30}"
  fi
  if [[ -z "${OPT_STARDIST_NMS_THRESH+x}" ]]; then
    OPT_STARDIST_NMS_THRESH="${STARDIST_NMS_THRESH:-}"
  fi
  if [[ -z "${OPT_STARDIST_MAX_DIM+x}" ]]; then
    OPT_STARDIST_MAX_DIM="${STARDIST_MAX_DIM:-6000}"
  fi
  if [[ -z "${OPT_STARDIST_SCALE+x}" ]]; then
    OPT_STARDIST_SCALE="${STARDIST_SCALE:-2.0}"
  fi
  if [[ -z "${OPT_STARDIST_N_TILES+x}" ]]; then
    OPT_STARDIST_N_TILES="${STARDIST_N_TILES:-4,4}"
  fi
  if [[ -z "${OPT_STARDIST_EXPAND_PIXELS+x}" ]]; then
    OPT_STARDIST_EXPAND_PIXELS="${STARDIST_EXPAND_PIXELS:-8}"
  fi
  if [[ -z "${OPT_STARDIST_MIN_UMI+x}" ]]; then
    OPT_STARDIST_MIN_UMI="${STARDIST_MIN_UMI:-10}"
  fi
  if [[ -z "${OPT_STARDIST_MIN_GENES+x}" ]]; then
    OPT_STARDIST_MIN_GENES="${STARDIST_MIN_GENES:-5}"
  fi
  if [[ -z "${OPT_STARDIST_SKIP_COUNTS+x}" ]]; then
    OPT_STARDIST_SKIP_COUNTS="${STARDIST_SKIP_COUNTS:-0}"
  fi
  if [[ -z "${OPT_STARDIST_NO_LABEL_OUTPUT+x}" ]]; then
    OPT_STARDIST_NO_LABEL_OUTPUT="${STARDIST_NO_LABEL_OUTPUT:-0}"
  fi
  if [[ -z "${OPT_STARDIST_TILED_INFERENCE+x}" ]]; then
    OPT_STARDIST_TILED_INFERENCE="${STARDIST_TILED_INFERENCE:-0}"
  fi
  if [[ -z "${OPT_STARDIST_TILE_SIZE+x}" ]]; then
    OPT_STARDIST_TILE_SIZE="${STARDIST_TILE_SIZE:-4096}"
  fi
  if [[ -z "${OPT_STARDIST_TILE_OVERLAP+x}" ]]; then
    OPT_STARDIST_TILE_OVERLAP="${STARDIST_TILE_OVERLAP:-256}"
  fi
  if [[ -z "${OPT_STARDIST_TILE_MERGE_OVERLAP+x}" ]]; then
    OPT_STARDIST_TILE_MERGE_OVERLAP="${STARDIST_TILE_MERGE_OVERLAP:-0.50}"
  fi
  if [[ -z "${OPT_STARDIST_FLUORESCENCE_CHANNEL+x}" ]]; then
    OPT_STARDIST_FLUORESCENCE_CHANNEL="${STARDIST_FLUORESCENCE_CHANNEL:-gray}"
  fi
  if [[ -z "${OPT_STARDIST_TILE_CACHE_DIR+x}" ]]; then
    OPT_STARDIST_TILE_CACHE_DIR="${STARDIST_TILE_CACHE_DIR:-}"
  fi
  if [[ -z "${OPT_STARDIST_EXCLUDE_RECTANGLES+x}" ]]; then
    OPT_STARDIST_EXCLUDE_RECTANGLES="${STARDIST_EXCLUDE_RECTANGLES:-}"
  fi
}

celatlas_post_binsegment_finalize_defaults() {
  local current_method="${1:-}"
  if [[ -n "${OPT_STARDIST_MODEL:-}" ]]; then
    return 0
  fi
  case "${current_method:-}" in
    image|ssDNA)
      OPT_STARDIST_MODEL="${STARDIST_MODEL:-2D_versatile_fluo}"
      ;;
    *)
      OPT_STARDIST_MODEL="${STARDIST_MODEL:-2D_versatile_he}"
      ;;
  esac
}

_celatlas_post_require_value() {
  local option="$1"
  local value="${2:-}"
  if [[ -z "$value" || "$value" == --* ]]; then
    echo "Error: ${option} requires a value"
    exit 1
  fi
}

celatlas_post_binsegment_parse_option() {
  CELATLAS_POST_ARG_SHIFT=0
  case "${1:-}" in
    --enable-cavity-filter)
      OPT_ENABLE_CAVITY_FILTER=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --disable-cavity-filter|--skip-cavity-filter)
      OPT_ENABLE_CAVITY_FILTER=0
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --cavity-bins)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_CAVITY_BINS="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --cavity-apply)
      OPT_CAVITY_APPLY=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --cavity-dry-run)
      OPT_CAVITY_DRY_RUN=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --cavity-skip-qc-images)
      OPT_CAVITY_SKIP_QC_IMAGES=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --cavity-mask-preset|--mask-preset)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_CAVITY_MASK_PRESET="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --enable-stardist-cell-segment|--enable-stardist)
      OPT_ENABLE_STARDIST_CELL_SEGMENT=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --disable-stardist-cell-segment|--skip-stardist|--disable-stardist)
      OPT_ENABLE_STARDIST_CELL_SEGMENT=0
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --stardist-python)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_PYTHON="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-labels)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_LABELS="$2"
      OPT_ENABLE_STARDIST_CELL_SEGMENT=1
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-image)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_IMAGE="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-tissue-bbox)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_TISSUE_BBOX="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-model)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_MODEL="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-model-dir)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_MODEL_DIR="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-prob-thresh)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_PROB_THRESH="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-nms-thresh)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_NMS_THRESH="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-max-dim)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_MAX_DIM="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-scale)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_SCALE="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-n-tiles)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_N_TILES="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-expand-pixels)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_EXPAND_PIXELS="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-min-umi)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_MIN_UMI="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-min-genes)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_MIN_GENES="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-skip-counts)
      OPT_STARDIST_SKIP_COUNTS=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --stardist-no-label-output)
      OPT_STARDIST_NO_LABEL_OUTPUT=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --stardist-tiled-inference)
      OPT_STARDIST_TILED_INFERENCE=1
      CELATLAS_POST_ARG_SHIFT=1
      return 0
      ;;
    --stardist-tile-size)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_TILE_SIZE="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-tile-overlap)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_TILE_OVERLAP="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-tile-merge-overlap)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_TILE_MERGE_OVERLAP="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-fluorescence-channel)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_FLUORESCENCE_CHANNEL="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-tile-cache-dir)
      _celatlas_post_require_value "$1" "${2:-}"
      OPT_STARDIST_TILE_CACHE_DIR="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
    --stardist-exclude-rectangle)
      _celatlas_post_require_value "$1" "${2:-}"
      if [[ -n "${OPT_STARDIST_EXCLUDE_RECTANGLES:-}" ]]; then
        OPT_STARDIST_EXCLUDE_RECTANGLES+=";"
      fi
      OPT_STARDIST_EXCLUDE_RECTANGLES+="$2"
      CELATLAS_POST_ARG_SHIFT=2
      return 0
      ;;
  esac
  return 1
}

celatlas_post_binsegment_help() {
  cat <<'EOF'
POST-BINSEGMENT OPTIONAL MODULES:
  --enable-cavity-filter        Run cavity/tissue mask filtering after 06.segment/01.binsegment
  --disable-cavity-filter       Disable cavity/tissue mask filtering
  --cavity-bins <bins>          Cavity-filtered bin list (default: 10,20,50,100,Raw)
  --cavity-apply                Replace 06.segment/01.binsegment/square_bin after backing it up
  --cavity-dry-run              Report cavity filtering counts only
  --cavity-skip-qc-images       Skip cavity-filter QC image generation
  --cavity-mask-preset <name>   Mask preset: default or adipose

  --enable-stardist-cell-segment
                                Run StarDist nucleus segmentation and barcode/count assignment
  --disable-stardist-cell-segment
                                Disable StarDist cell segmentation
  --stardist-python <path>      Python executable with stardist/tensorflow installed
  --stardist-labels <path>      Existing StarDist label TIFF; skips inference
  --stardist-image <path>       DAPI/fluorescence input or matching registered preview image
  --stardist-tissue-bbox <path> Label canvas bbox (x0,y0,width,height) for barcode assignment
  --stardist-model <name>       StarDist model name (default: HE -> 2D_versatile_he; image/ssDNA -> 2D_versatile_fluo)
  --stardist-model-dir <path>   Local StarDist model directory
  --stardist-prob-thresh <num>  StarDist probability threshold (default: 0.30)
  --stardist-nms-thresh <num>   Optional StarDist NMS threshold
  --stardist-max-dim <num>      Inference image long edge (default: 6000)
  --stardist-scale <num>        StarDist internal scale (default: 2.0)
  --stardist-n-tiles <grid>     Tile grid, e.g. 4,4
  --stardist-expand-pixels <n>  Label expansion radius before barcode assignment
  --stardist-min-umi <n>        Minimum UMI per retained cell
  --stardist-min-genes <n>      Minimum genes per retained cell
  --stardist-skip-counts        Skip count_detail aggregation
  --stardist-no-label-output    Do not write slide-sized expanded/QC label TIFFs
  --stardist-tiled-inference    Run full-resolution overlapping tile inference
  --stardist-tile-size <n>      Disk tile size (default: 4096)
  --stardist-tile-overlap <n>   Tile context halo (default: 256)
  --stardist-tile-merge-overlap <num>
                                Duplicate merge fraction (default: 0.50)
  --stardist-fluorescence-channel <name>
                                gray, red, green, or blue (default: gray)
  --stardist-tile-cache-dir <path>
                                Resumable tile cache directory
  --stardist-exclude-rectangle <x0,y0,x1,y1>
                                Remove annotation/scale-bar region; repeatable

EOF
}

_celatlas_post_run_command() {
  if declare -f run_command >/dev/null 2>&1; then
    run_command "$@"
  else
    echo "Running:" "$@"
    "$@"
  fi
}

celatlas_run_post_binsegment_optional_modules() {
  celatlas_post_binsegment_init_defaults

  if [[ "${OPT_ENABLE_CAVITY_FILTER}" != "1" && "${OPT_ENABLE_STARDIST_CELL_SEGMENT}" != "1" ]]; then
    return 0
  fi

  if [[ -z "${runtime_mask_dir:-}" ]]; then
    runtime_mask_dir="${sampledir}/06.segment/mask"
  fi

  local required_name
  for required_name in sampledir sample runtime_mask_dir method; do
    if [[ -z "${!required_name:-}" ]]; then
      echo "ERROR: ${required_name} is required before running post-binSegment modules"
      exit 1
    fi
  done

  celatlas_post_binsegment_finalize_defaults "$method"

  if [[ "$OPT_ENABLE_CAVITY_FILTER" == "1" ]]; then
    echo ""
    echo "======================================"
    echo "Step 06.segment/01: Cavity Filter"
    echo "======================================"

    local cavity_args=(
      rna cavityFilter
      --sampledir "${sampledir}"
      --bin-segment-dir "${sampledir}/06.segment/01.binsegment"
      --sample "${sample}"
      --bins "${OPT_CAVITY_BINS}"
      --mask-preset "${OPT_CAVITY_MASK_PRESET}"
    )
    if [[ "$OPT_CAVITY_APPLY" == "1" ]]; then
      cavity_args+=(--apply)
    fi
    if [[ "$OPT_CAVITY_DRY_RUN" == "1" ]]; then
      cavity_args+=(--dry-run)
    fi
    if [[ "$OPT_CAVITY_SKIP_QC_IMAGES" == "1" ]]; then
      cavity_args+=(--skip-qc-images)
    fi
    _celatlas_post_run_command celatlas_spatial "${cavity_args[@]}"
  fi

  if [[ "$OPT_ENABLE_STARDIST_CELL_SEGMENT" == "1" ]]; then
    echo ""
    echo "======================================"
    echo "Step 06.segment/02: StarDist Cell Segment"
    echo "======================================"
    if [[ "$method" != "HE" ]]; then
      echo "WARNING: StarDist H&E segmentation is intended for method=HE; current method=${method}."
    fi

    local stardist_args=(
      rna stardistCellSegment
      --sampledir "${sampledir}"
      --bin-segment-dir "${sampledir}/06.segment/01.binsegment"
      --outdir "${sampledir}/06.segment/02.cellsegment"
      --sample "${sample}"
      --genomeDir "${genomeDir}"
      --runtime-input-dir "${runtime_mask_dir}"
      --model "${OPT_STARDIST_MODEL}"
      --prob-thresh "${OPT_STARDIST_PROB_THRESH}"
      --max-dim "${OPT_STARDIST_MAX_DIM}"
      --scale "${OPT_STARDIST_SCALE}"
      --n-tiles "${OPT_STARDIST_N_TILES}"
      --expand-pixels "${OPT_STARDIST_EXPAND_PIXELS}"
      --min-umi "${OPT_STARDIST_MIN_UMI}"
      --min-genes "${OPT_STARDIST_MIN_GENES}"
    )
    if [[ -n "$OPT_STARDIST_PYTHON" ]]; then
      stardist_args+=(--stardist-python "$OPT_STARDIST_PYTHON")
    fi
    if [[ -n "$OPT_STARDIST_LABELS" ]]; then
      stardist_args+=(--labels "$OPT_STARDIST_LABELS" --skip-inference)
    fi
    if [[ -n "$OPT_STARDIST_IMAGE" ]]; then
      stardist_args+=(--image "$OPT_STARDIST_IMAGE")
    fi
    if [[ -n "$OPT_STARDIST_TISSUE_BBOX" ]]; then
      stardist_args+=(--tissue-bbox "$OPT_STARDIST_TISSUE_BBOX")
    fi
    if [[ -n "$OPT_STARDIST_MODEL_DIR" ]]; then
      stardist_args+=(--model-dir "$OPT_STARDIST_MODEL_DIR")
    fi
    if [[ -n "$OPT_STARDIST_NMS_THRESH" ]]; then
      stardist_args+=(--nms-thresh "$OPT_STARDIST_NMS_THRESH")
    fi
    if [[ "$OPT_STARDIST_SKIP_COUNTS" == "1" ]]; then
      stardist_args+=(--skip-counts)
    fi
    if [[ "$OPT_STARDIST_NO_LABEL_OUTPUT" == "1" ]]; then
      stardist_args+=(--no-label-output)
    fi
    if [[ "$OPT_STARDIST_TILED_INFERENCE" == "1" ]]; then
      stardist_args+=(--tiled-inference)
    fi
    if [[ -n "$OPT_STARDIST_TILE_SIZE" ]]; then
      stardist_args+=(--tile-size "$OPT_STARDIST_TILE_SIZE")
    fi
    if [[ -n "$OPT_STARDIST_TILE_OVERLAP" ]]; then
      stardist_args+=(--tile-overlap "$OPT_STARDIST_TILE_OVERLAP")
    fi
    if [[ -n "$OPT_STARDIST_TILE_MERGE_OVERLAP" ]]; then
      stardist_args+=(--tile-merge-overlap "$OPT_STARDIST_TILE_MERGE_OVERLAP")
    fi
    if [[ -n "$OPT_STARDIST_FLUORESCENCE_CHANNEL" ]]; then
      stardist_args+=(--fluorescence-channel "$OPT_STARDIST_FLUORESCENCE_CHANNEL")
    fi
    if [[ -n "$OPT_STARDIST_TILE_CACHE_DIR" ]]; then
      stardist_args+=(--tile-cache-dir "$OPT_STARDIST_TILE_CACHE_DIR")
    fi
    if [[ -n "${OPT_STARDIST_EXCLUDE_RECTANGLES:-}" ]]; then
      local old_ifs="$IFS"
      IFS=';'
      read -r -a exclude_rectangles <<< "$OPT_STARDIST_EXCLUDE_RECTANGLES"
      IFS="$old_ifs"
      local rectangle
      for rectangle in "${exclude_rectangles[@]}"; do
        [[ -n "$rectangle" ]] && stardist_args+=(--exclude-rectangle "$rectangle")
      done
    fi
    _celatlas_post_run_command celatlas_spatial "${stardist_args[@]}"
  fi
}
