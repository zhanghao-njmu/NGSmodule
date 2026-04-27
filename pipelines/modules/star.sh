#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/star.sh — RNA-seq alignment with STAR.
#
# Functions:
#   module_star_align --index <dir> --out-prefix <p>
#                     --in1 <R1> [--in2 <R2>]
#                     [--threads N] [--sam-type "BAM SortedByCoordinate"]
#                     [--read-cmd zcat]
#
# Output convention: ${out-prefix}Aligned.sortedByCoord.out.bam (STAR's
# default naming when --outSAMtype BAM SortedByCoordinate). The pipeline
# is responsible for renaming/symlinking to ${sample}.STAR.sorted.bam.
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_star_align() {
  local index="" out_prefix="" in1="" in2="" extra=""
  local threads="${threads:-4}"
  local sam_type="BAM SortedByCoordinate"
  local read_cmd="zcat"
  while (( $# > 0 )); do
    case "$1" in
      --index)       index="$2"; shift 2 ;;
      --out-prefix)  out_prefix="$2"; shift 2 ;;
      --in1)         in1="$2"; shift 2 ;;
      --in2)         in2="$2"; shift 2 ;;
      --threads)     threads="$2"; shift 2 ;;
      --sam-type)    sam_type="$2"; shift 2 ;;
      --read-cmd)    read_cmd="$2"; shift 2 ;;
      --extra)       extra="$2"; shift 2 ;;
      *) echo "module_star_align: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${index:?--index required}"; : "${out_prefix:?--out-prefix required}"
  : "${in1:?--in1 required}"

  local reads_arg=("$in1")
  [[ -n "$in2" ]] && reads_arg+=("$in2")

  # shellcheck disable=SC2086  # sam_type and extra intentionally split
  module_run "STAR align" \
    STAR \
      --runThreadN "$threads" \
      --genomeDir "$index" \
      --readFilesIn "${reads_arg[@]}" \
      --readFilesCommand "$read_cmd" \
      --outFileNamePrefix "$out_prefix" \
      --outSAMtype $sam_type \
      $extra
}
