#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/bwa_mem2.sh — DNA short-read alignment with BWA-MEM2.
#
# Outputs unsorted SAM on stdout; pair with module_samtools_sort to get
# a coordinate-sorted BAM (this is the conventional pipeline pattern and
# matches the bwa-mem2 docs).
#
#   module_bwa_mem2_align --index <prefix> --in1 <R1> [--in2 <R2>]
#                         --rg-sample <sm> [--rg-id <id>] [--rg-lib <lib>]
#                         [--threads N]
#                         | module_samtools_sort --out file.bam
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_bwa_mem2_align() {
  local index="" in1="" in2="" rg_id="" rg_sm="" rg_lib=""
  local threads="${threads:-4}"
  while (( $# > 0 )); do
    case "$1" in
      --index)      index="$2"; shift 2 ;;
      --in1)        in1="$2"; shift 2 ;;
      --in2)        in2="$2"; shift 2 ;;
      --rg-id)      rg_id="$2"; shift 2 ;;
      --rg-sample)  rg_sm="$2"; shift 2 ;;
      --rg-lib)     rg_lib="$2"; shift 2 ;;
      --threads)    threads="$2"; shift 2 ;;
      *) echo "module_bwa_mem2_align: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${index:?--index required}"; : "${in1:?--in1 required}"
  : "${rg_sm:?--rg-sample required}"
  rg_id="${rg_id:-$rg_sm}"
  rg_lib="${rg_lib:-$rg_sm}"

  local rg="@RG\tID:${rg_id}\tSM:${rg_sm}\tLB:${rg_lib}\tPL:ILLUMINA"
  local reads=("$in1")
  [[ -n "$in2" ]] && reads+=("$in2")

  module_run "bwa-mem2 align" \
    bwa-mem2 mem -t "$threads" -R "$rg" "$index" "${reads[@]}"
}
