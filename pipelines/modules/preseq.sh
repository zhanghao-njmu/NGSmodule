#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/preseq.sh — Library complexity estimation.
#
# Functions:
#   module_preseq_lc_extrap --bam <bam> --out <txt> [--extra "..."]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_preseq_lc_extrap() {
  local bam="" out="" extra=""
  while (( $# > 0 )); do
    case "$1" in
      --bam)   bam="$2"; shift 2 ;;
      --out)   out="$2"; shift 2 ;;
      --extra) extra="$2"; shift 2 ;;
      *) echo "module_preseq_lc_extrap: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${out:?--out required}"

  # shellcheck disable=SC2086  # extra is intentionally word-split
  module_run "preseq lc_extrap" \
    preseq lc_extrap -B "$bam" -o "$out" $extra
}
