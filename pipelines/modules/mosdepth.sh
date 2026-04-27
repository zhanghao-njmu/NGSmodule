#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/mosdepth.sh — Fast BAM/CRAM depth calculation.
#
# Functions:
#   module_mosdepth_run --bam <bam> --prefix <p>
#                       [--threads N] [--no-per-base] [--fast-mode]
#                       [--by <bed>] [--extra "..."]
#
# Default flags mirror the legacy postAlignmentQC: -n --fast-mode.
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_mosdepth_run() {
  local bam="" prefix="" by="" extra=""
  local threads="${threads:-4}"
  local no_per_base=1
  local fast_mode=1
  while (( $# > 0 )); do
    case "$1" in
      --bam)         bam="$2"; shift 2 ;;
      --prefix)      prefix="$2"; shift 2 ;;
      --threads)     threads="$2"; shift 2 ;;
      --by)          by="$2"; shift 2 ;;
      --extra)       extra="$2"; shift 2 ;;
      --no-per-base) no_per_base=1; shift ;;
      --per-base)    no_per_base=0; shift ;;
      --fast-mode)   fast_mode=1; shift ;;
      --no-fast-mode) fast_mode=0; shift ;;
      *) echo "module_mosdepth_run: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${prefix:?--prefix required}"

  local args=(-t "$threads")
  (( no_per_base )) && args+=(-n)
  (( fast_mode ))   && args+=(--fast-mode)
  [[ -n "$by" ]]    && args+=(--by "$by")

  # shellcheck disable=SC2086  # extra is intentionally word-split
  module_run "mosdepth" \
    mosdepth "${args[@]}" $extra "$prefix" "$bam"
}
