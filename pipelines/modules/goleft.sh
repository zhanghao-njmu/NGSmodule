#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/goleft.sh — Coverage analysis from BAM index.
#
# Functions:
#   module_goleft_indexcov --bam <bam> --dir <out_dir> [--extra "..."]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_goleft_indexcov() {
  local bam="" dir="" extra=""
  while (( $# > 0 )); do
    case "$1" in
      --bam)   bam="$2"; shift 2 ;;
      --dir)   dir="$2"; shift 2 ;;
      --extra) extra="$2"; shift 2 ;;
      *) echo "module_goleft_indexcov: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${dir:?--dir required}"
  mkdir -p "$dir"

  # shellcheck disable=SC2086  # extra is intentionally word-split
  module_run "goleft indexcov" \
    goleft indexcov --directory "$dir" $extra "$bam"
}
