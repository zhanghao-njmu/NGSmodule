#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/samtools.sh — sort, index, flagstat.
#
# Functions:
#   module_samtools_sort     --in <sam|bam|->  --out <bam>  [--threads N]
#   module_samtools_index    --bam <bam>
#   module_samtools_flagstat --bam <bam>  --out <txt>
#
# `--in -` reads from stdin (so you can pipe an aligner directly into sort).
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_samtools_sort() {
  local in="" out="" threads="${threads:-4}"
  while (( $# > 0 )); do
    case "$1" in
      --in)      in="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *) echo "module_samtools_sort: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"

  module_run "samtools sort" \
    samtools sort -@ "$threads" -o "$out" "$in"
}

module_samtools_index() {
  local bam=""
  while (( $# > 0 )); do
    case "$1" in
      --bam) bam="$2"; shift 2 ;;
      *) echo "module_samtools_index: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"
  module_run "samtools index" samtools index "$bam"
}

module_samtools_flagstat() {
  local bam="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --bam) bam="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      *) echo "module_samtools_flagstat: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${out:?--out required}"

  if declare -f is_dry_run >/dev/null && is_dry_run; then
    printf '::dry-run::samtools flagstat %s > %s\n' "$bam" "$out"
    NGS_MODULE_CURRENT_TOOL="samtools" module_run "samtools flagstat" true
    return 0
  fi
  NGS_MODULE_CURRENT_TOOL="samtools" \
    module_run "samtools flagstat" \
    bash -c "samtools flagstat '$bam' > '$out'"
}
