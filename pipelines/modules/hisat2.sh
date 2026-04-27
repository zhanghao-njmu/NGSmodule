#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/hisat2.sh — splice-aware RNA aligner.
#
#   module_hisat2_align --index <prefix> --in1 <R1> [--in2 <R2>]
#                       [--threads N] [--strandness {RF,FR,unstranded}]
#                       --out-sam <path>
#
# Writes unsorted SAM (pair with samtools sort to get sorted BAM).
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_hisat2_align() {
  local index="" in1="" in2="" out_sam=""
  local threads="${threads:-4}" strand=""
  while (( $# > 0 )); do
    case "$1" in
      --index)      index="$2"; shift 2 ;;
      --in1)        in1="$2"; shift 2 ;;
      --in2)        in2="$2"; shift 2 ;;
      --out-sam)    out_sam="$2"; shift 2 ;;
      --threads)    threads="$2"; shift 2 ;;
      --strandness) strand="$2"; shift 2 ;;
      *) echo "module_hisat2_align: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${index:?--index required}"; : "${in1:?--in1 required}"
  : "${out_sam:?--out-sam required}"

  local strand_arg=()
  case "${strand,,}" in
    rf) strand_arg=(--rna-strandness RF) ;;
    fr) strand_arg=(--rna-strandness FR) ;;
    *)  : ;;
  esac

  if [[ -n "$in2" ]]; then
    module_run "hisat2 PE" \
      hisat2 -p "$threads" -x "$index" \
        -1 "$in1" -2 "$in2" -S "$out_sam" "${strand_arg[@]}"
  else
    module_run "hisat2 SE" \
      hisat2 -p "$threads" -x "$index" \
        -U "$in1" -S "$out_sam" "${strand_arg[@]}"
  fi
}
