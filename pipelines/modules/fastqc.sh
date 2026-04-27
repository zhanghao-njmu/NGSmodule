#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/fastqc.sh
#
#   module_fastqc_run --out-dir <dir> [--threads N] <fastq>...
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_fastqc_run() {
  local out_dir="" threads="${threads:-4}"
  local files=()
  while (( $# > 0 )); do
    case "$1" in
      --out-dir) out_dir="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *)         files+=("$1"); shift ;;
    esac
  done
  : "${out_dir:?--out-dir required}"
  if (( ${#files[@]} == 0 )); then
    echo "module_fastqc_run: at least one fastq required" >&2; return 2
  fi
  mkdir -p "$out_dir"

  module_run "fastqc" \
    fastqc -o "$out_dir" -t "$threads" --quiet "${files[@]}"
}
