#!/usr/bin/env bash
set -euo pipefail

# Unified Celatlas v1.8 entrypoint.
# This wrapper only selects and validates the target pipeline script. It keeps
# Celatlas.sh, Celatlas_FFPE.sh, Celatlas_SN.sh, and Celatlas_reanalysis.sh as
# the source of truth for the real pipeline steps.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

show_help() {
  cat <<'EOF'
=====================================================================
Celatlas v1.8 Unified Runner
=====================================================================

USAGE:
  bash Celatlas_run.sh <workflow> <pipeline arguments...>
  bash Celatlas_run.sh --workflow <workflow> <pipeline arguments...>

WORKFLOWS:
  CONFIG, config
      Read a per-job config file and translate it to FF/FFPE/SN/reanalysis.
      This is the preferred production integration mode.

  FF, fresh, standard
      Run Celatlas.sh for the standard fresh/frozen spatial workflow.

  FFPE, probe, panel
      Run Celatlas_FFPE.sh for FFPE / targeted-panel / H&E workflows.

  SN, single-nucleus, snrna
      Run Celatlas_SN.sh for single-nucleus RNA workflows.

  reanalysis, re
      Run Celatlas_reanalysis.sh from existing 05.count output, typically
      after manual registration or mask adjustment.

WRAPPER OPTIONS:
  --workflow, --pipeline, --type <workflow>
      Select FF, FFPE, SN, or reanalysis without using the first positional token.

  --dry-run
      Print the selected target command without running it.

  --help, -h
      Show this wrapper help. If placed after a workflow, it is passed to the
      selected target script, for example: bash Celatlas_run.sh FFPE --help

EXAMPLES:
  # Production config mode
  bash Celatlas_run.sh CONFIG configs/jobs/example.strna.env --dry-run

  # Standard fresh/frozen run
  bash Celatlas_run.sh FF test STHD110017_A1 hd BBV4 Mus_musculus gene_expr strna

  # FFPE H&E full run with cavity filtering and StarDist cell segmentation
  bash Celatlas_run.sh FFPE test SX000258_B1 SX_test BBV2.4 Homo_sapiens HE strna \
    --thread 32 --enable-cavity-filter --cavity-apply --enable-stardist-cell-segment

  # Single-nucleus run
  bash Celatlas_run.sh SN test SX000256_A1 SN_test BBV2.4 Sorghum_bicolorv_geneid gene_expr strna \
    --thread 32 --bin 50

  # Reanalysis after manual registration/mask adjustment
  bash Celatlas_run.sh reanalysis SX000258_B1 SX_test BBV2.4 Homo_sapiens HE strna \
    --thread 32 --enable-cavity-filter --cavity-apply --enable-stardist-cell-segment

  # Dry run to check which script and arguments will be used
  bash Celatlas_run.sh --workflow FFPE test SX000258_B1 SX_test BBV2.4 Homo_sapiens HE strna --dry-run

NOTES:
  - Use --workflow/--pipeline/--type for the wrapper selector.
  - Keep --mode for the underlying Celatlas assay mode, usually strna or scrna.
  - This wrapper requires an explicit workflow to avoid accidental use of the
    wrong production script.
EOF
}

die() {
  echo "ERROR: $*" >&2
  echo "Run 'bash ${BASH_SOURCE[0]} --help' for usage information." >&2
  exit 1
}

warn() {
  echo "WARNING: $*" >&2
}

lowercase() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

normalize_workflow() {
  local raw
  raw="$(lowercase "${1:-}")"
  case "$raw" in
    config|conf|cfg|job|env)
      printf 'CONFIG'
      ;;
    ff|fresh|fresh-frozen|fresh_frozen|standard|std|normal)
      printf 'FF'
      ;;
    ffpe|probe|panel|targeted|targeted-panel|targeted_panel|wtpanel)
      printf 'FFPE'
      ;;
    sn|single-nucleus|single_nucleus|snrna|sn-rna|sn_rna)
      printf 'SN'
      ;;
    reanalysis|reanalyze|reanalyse|re|manual|manual-registration|manual_registration)
      printf 'reanalysis'
      ;;
    *)
      return 1
      ;;
  esac
}

target_for_workflow() {
  case "$1" in
    CONFIG)
      printf '%s/scripts/celatlas_config_run.sh' "$SCRIPT_DIR"
      ;;
    FF)
      printf '%s/Celatlas.sh' "$SCRIPT_DIR"
      ;;
    FFPE)
      printf '%s/Celatlas_FFPE.sh' "$SCRIPT_DIR"
      ;;
    SN)
      printf '%s/Celatlas_SN.sh' "$SCRIPT_DIR"
      ;;
    reanalysis)
      printf '%s/Celatlas_reanalysis.sh' "$SCRIPT_DIR"
      ;;
    *)
      return 1
      ;;
  esac
}

is_value_option() {
  case "$1" in
    --workflow|--pipeline|--type|\
    --sample|-s|--chip|-c|--casno|--chemistry|--species|--method|-m|--mode|\
    --thread|-t|--bin|-b|--insertR2|--cell_num|--pixelSize|\
    --bbv4-strna-barcode-mismatch|--bbv4-strna-mismatch|\
    --star-match-min|--star-match-ratio|--star-score-ratio|\
    --reference_dir|--mask_dir|--image_dir|--fastq_dir|--fastq_name|\
    --star-multimap|--cell-cluster-resolution|--cell-n-neighbors|--cell-n-pcs|\
    --workspace|--genomeDir|--sampledir|--src_dir|\
    --cavity-bins|--cavity-mask-preset|\
    --stardist-python|--stardist-labels|--stardist-image|--stardist-tissue-bbox|--stardist-model|--stardist-model-dir|\
    --stardist-prob-thresh|--stardist-nms-thresh|--stardist-max-dim|\
    --stardist-scale|--stardist-n-tiles|--stardist-expand-pixels|\
    --stardist-min-umi|--stardist-min-genes|\
    --stardist-tile-size|--stardist-tile-overlap|--stardist-tile-merge-overlap|\
    --stardist-fluorescence-channel|--stardist-tile-cache-dir|--stardist-exclude-rectangle)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

is_flag_option() {
  case "$1" in
    --enable-cavity-filter|--disable-cavity-filter|--skip-cavity-filter|\
    --cavity-apply|--cavity-dry-run|--cavity-skip-qc-images|\
    --enable-stardist-cell-segment|--enable-stardist|\
    --disable-stardist-cell-segment|--skip-stardist|--disable-stardist|\
    --stardist-skip-counts|--stardist-no-label-output|--stardist-tiled-inference|\
    --resume-existing|\
    --resolve-multigene-umi|\
    --no-clean|--skip-binsegment|--skip-analysis|--skip-report|\
    --gene-mask-filter|\
    --debug)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

parse_target_summary() {
  local workflow="$1"
  shift

  SUMMARY_SAMPLE=""
  SUMMARY_CHIP=""
  SUMMARY_CASNO=""
  SUMMARY_CHEMISTRY=""
  SUMMARY_SPECIES=""
  SUMMARY_METHOD=""
  SUMMARY_MODE=""
  SUMMARY_THREAD=""
  SUMMARY_BIN=""
  SUMMARY_SAMPLEDIR=""
  SUMMARY_NO_CLEAN=0
  SUMMARY_SKIP_BINSEGMENT=0
  SUMMARY_SKIP_ANALYSIS=0
  SUMMARY_SKIP_REPORT=0
  SUMMARY_CAVITY=0
  SUMMARY_CAVITY_APPLY=0
  SUMMARY_STARDIST=0
  SUMMARY_GENE_MASK_FILTER=0

  local positional=()
  local arg next
  while (($# > 0)); do
    arg="$1"
    case "$arg" in
      --sample|-s)
        SUMMARY_SAMPLE="${2:-}"
        shift 2
        ;;
      --chip|-c)
        SUMMARY_CHIP="${2:-}"
        shift 2
        ;;
      --casno)
        SUMMARY_CASNO="${2:-}"
        shift 2
        ;;
      --chemistry)
        SUMMARY_CHEMISTRY="${2:-}"
        shift 2
        ;;
      --species)
        SUMMARY_SPECIES="${2:-}"
        shift 2
        ;;
      --method|-m)
        SUMMARY_METHOD="${2:-}"
        shift 2
        ;;
      --mode)
        SUMMARY_MODE="${2:-}"
        shift 2
        ;;
      --thread|-t)
        SUMMARY_THREAD="${2:-}"
        shift 2
        ;;
      --bin|-b)
        SUMMARY_BIN="${2:-}"
        shift 2
        ;;
      --sampledir)
        SUMMARY_SAMPLEDIR="${2:-}"
        shift 2
        ;;
      --no-clean)
        SUMMARY_NO_CLEAN=1
        shift
        ;;
      --skip-binsegment)
        SUMMARY_SKIP_BINSEGMENT=1
        shift
        ;;
      --skip-analysis)
        SUMMARY_SKIP_ANALYSIS=1
        shift
        ;;
      --skip-report)
        SUMMARY_SKIP_REPORT=1
        shift
        ;;
      --enable-cavity-filter)
        SUMMARY_CAVITY=1
        shift
        ;;
      --cavity-apply)
        SUMMARY_CAVITY_APPLY=1
        shift
        ;;
      --enable-stardist-cell-segment|--enable-stardist)
        SUMMARY_STARDIST=1
        shift
        ;;
      --gene-mask-filter)
        SUMMARY_GENE_MASK_FILTER=1
        shift
        ;;
      --)
        shift
        while (($# > 0)); do
          positional+=("$1")
          shift
        done
        ;;
      -*)
        if is_value_option "$arg"; then
          [[ $# -ge 2 ]] || die "Option '$arg' requires a value"
          shift 2
        elif is_flag_option "$arg"; then
          shift
        else
          # Unknown options are intentionally passed through. If they take a
          # value, the target script will validate them.
          shift
        fi
        ;;
      *)
        positional+=("$arg")
        shift
        ;;
    esac
  done

  local n=${#positional[@]}
  if [[ "$workflow" == "FF" ]]; then
    if (( n > 0 && n != 6 && n != 7 && n != 11 && n != 12 )); then
      die "FF/Celatlas.sh positional arguments must be 6, 7, 11, or 12; got $n"
    fi
  elif (( n > 0 && n != 6 && n != 7 )); then
    die "${workflow}/target positional arguments must be 6 or 7; got $n"
  fi

  if (( n >= 7 )); then
    SUMMARY_SAMPLE="${SUMMARY_SAMPLE:-${positional[0]}}"
    SUMMARY_CHIP="${SUMMARY_CHIP:-${positional[1]}}"
    SUMMARY_CASNO="${SUMMARY_CASNO:-${positional[2]}}"
    SUMMARY_CHEMISTRY="${SUMMARY_CHEMISTRY:-${positional[3]}}"
    SUMMARY_SPECIES="${SUMMARY_SPECIES:-${positional[4]}}"
    SUMMARY_METHOD="${SUMMARY_METHOD:-${positional[5]}}"
    SUMMARY_MODE="${SUMMARY_MODE:-${positional[6]}}"
  elif (( n == 6 )); then
    SUMMARY_CHIP="${SUMMARY_CHIP:-${positional[0]}}"
    SUMMARY_CASNO="${SUMMARY_CASNO:-${positional[1]}}"
    SUMMARY_CHEMISTRY="${SUMMARY_CHEMISTRY:-${positional[2]}}"
    SUMMARY_SPECIES="${SUMMARY_SPECIES:-${positional[3]}}"
    SUMMARY_METHOD="${SUMMARY_METHOD:-${positional[4]}}"
    SUMMARY_MODE="${SUMMARY_MODE:-${positional[5]}}"
  fi

  if [[ -n "$SUMMARY_METHOD" ]]; then
    case "$SUMMARY_METHOD" in
      image|ssDNA|gene_expr|HE) ;;
      *) die "Invalid method '$SUMMARY_METHOD'. Expected image, ssDNA, gene_expr, or HE." ;;
    esac
  fi

  if [[ -n "$SUMMARY_MODE" ]]; then
    case "$SUMMARY_MODE" in
      strna|scrna) ;;
      *) die "Invalid assay mode '$SUMMARY_MODE'. Expected strna or scrna." ;;
    esac
  fi

  if [[ "$workflow" == "FFPE" && -n "$SUMMARY_METHOD" && "$SUMMARY_METHOD" != "HE" ]]; then
    warn "Workflow FFPE is usually paired with method HE; current method is '$SUMMARY_METHOD'."
  fi

  if [[ "$workflow" == "SN" && -n "$SUMMARY_MODE" && "$SUMMARY_MODE" != "strna" && "$SUMMARY_MODE" != "scrna" ]]; then
    warn "Workflow SN usually uses mode strna or scrna; current mode is '$SUMMARY_MODE'."
  fi

  if [[ "$workflow" == "reanalysis" && "$SUMMARY_NO_CLEAN" != "1" ]]; then
    warn "reanalysis will clean downstream 06.segment/07.outs unless you add --no-clean or skip flags."
  fi
}

print_summary() {
  local workflow="$1"
  local target_script="$2"

  echo "==============================================="
  echo "Celatlas Unified Runner"
  echo "==============================================="
  echo "Workflow:          $workflow"
  echo "Target Script:     $target_script"
  echo "Sample Name:       ${SUMMARY_SAMPLE:-N/A}"
  echo "Chip Number:       ${SUMMARY_CHIP:-N/A}"
  echo "Case Number:       ${SUMMARY_CASNO:-N/A}"
  echo "Chemistry:         ${SUMMARY_CHEMISTRY:-N/A}"
  echo "Species:           ${SUMMARY_SPECIES:-N/A}"
  echo "Method:            ${SUMMARY_METHOD:-N/A}"
  echo "Mode:              ${SUMMARY_MODE:-N/A}"
  echo "Thread:            ${SUMMARY_THREAD:-target default}"
  echo "Bin:               ${SUMMARY_BIN:-target default}"
  if [[ -n "$SUMMARY_SAMPLEDIR" ]]; then
    echo "Sample Dir:        $SUMMARY_SAMPLEDIR"
  fi
  echo "Cavity Filter:     $SUMMARY_CAVITY"
  echo "Cavity Apply:      $SUMMARY_CAVITY_APPLY"
  echo "StarDist Segment:  $SUMMARY_STARDIST"
  echo "Gene Mask Filter:  $SUMMARY_GENE_MASK_FILTER"
  if [[ "$workflow" == "reanalysis" ]]; then
    echo "No Clean:          $SUMMARY_NO_CLEAN"
    echo "Skip BinSegment:   $SUMMARY_SKIP_BINSEGMENT"
    echo "Skip Analysis:     $SUMMARY_SKIP_ANALYSIS"
    echo "Skip Report:       $SUMMARY_SKIP_REPORT"
  fi
  echo "==============================================="
}

quote_command() {
  printf '%q ' "$@"
  printf '\n'
}

workflow_raw=""
workflow=""
dry_run=0
target_args=()

while (($# > 0)); do
  case "$1" in
    --help|-h)
      if [[ -z "$workflow_raw" ]]; then
        show_help
        exit 0
      fi
      target_args+=("$1")
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    --workflow|--pipeline|--type)
      [[ $# -ge 2 ]] || die "Option '$1' requires a workflow value"
      workflow_raw="$2"
      shift 2
      ;;
    --)
      shift
      while (($# > 0)); do
        target_args+=("$1")
        shift
      done
      ;;
    *)
      if [[ -z "$workflow_raw" && ${#target_args[@]} -eq 0 ]]; then
        if normalize_workflow "$1" >/dev/null 2>&1; then
          workflow_raw="$1"
          shift
        else
          die "First argument must be a workflow: CONFIG, FF, FFPE, or reanalysis. Got '$1'"
        fi
      else
        target_args+=("$1")
        shift
      fi
      ;;
  esac
done

[[ -n "$workflow_raw" ]] || die "Missing workflow. Use CONFIG, FF, FFPE, SN, or reanalysis."
workflow="$(normalize_workflow "$workflow_raw")" || die "Unknown workflow '$workflow_raw'. Use CONFIG, FF, FFPE, SN, or reanalysis."
target_script="$(target_for_workflow "$workflow")" || die "No target script for workflow '$workflow'"

[[ -f "$target_script" ]] || die "Target script not found: $target_script"
[[ -r "$target_script" ]] || die "Target script is not readable: $target_script"

if [[ "$workflow" == "CONFIG" ]]; then
  echo "==============================================="
  echo "Celatlas Unified Runner"
  echo "==============================================="
  echo "Workflow:          $workflow"
  echo "Target Script:     $target_script"
  echo "Config Args:       ${target_args[*]:-N/A}"
  echo "==============================================="
else
  parse_target_summary "$workflow" "${target_args[@]}"
  print_summary "$workflow" "$target_script"
fi

cmd=(bash "$target_script" "${target_args[@]}")
if [[ "$workflow" == "CONFIG" && "$dry_run" == "1" ]]; then
  cmd+=(--dry-run)
fi
echo "Command:"
quote_command "${cmd[@]}"

if [[ "$dry_run" == "1" ]]; then
  if [[ "$workflow" == "CONFIG" ]]; then
    exec "${cmd[@]}"
  fi
  echo "[DRY-RUN] Command not executed."
  exit 0
fi

exec "${cmd[@]}"
