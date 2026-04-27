#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/bismark.sh — Bismark bisulfite alignment & extraction.
#
# Functions:
#   module_bismark_methylation_extractor
#       --bam <bam> --out-dir <dir> --genome-folder <dir> [--threads N]
#       [--paired] [--single]
#       [--extra "..."]
#
#   module_bismark2report --dir <dir>
#                         --alignment-report <file>
#                         [--splitting-report <file>] [--mbias-report <file>]
#                         [--dedup-report <file>] [--nucleotide-report <file>]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_bismark_methylation_extractor() {
  local bam="" out_dir="" genome_folder="" extra=""
  local threads="${threads:-4}"
  local layout="auto"   # auto | paired | single
  while (( $# > 0 )); do
    case "$1" in
      --bam)            bam="$2"; shift 2 ;;
      --out-dir)        out_dir="$2"; shift 2 ;;
      --genome-folder)  genome_folder="$2"; shift 2 ;;
      --threads)        threads="$2"; shift 2 ;;
      --paired)         layout="paired"; shift ;;
      --single)         layout="single"; shift ;;
      --extra)          extra="$2"; shift 2 ;;
      *) echo "module_bismark_methylation_extractor: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${out_dir:?--out-dir required}"
  : "${genome_folder:?--genome-folder required}"

  # bismark_methylation_extractor uses --multicore N (parallel chrom workers,
  # each spawning ~3 child processes); leave headroom by passing N/3 minimum 1.
  local mc=$(( threads / 3 ))
  (( mc < 1 )) && mc=1

  local args=(
    --multicore "$mc"
    --gzip
    --comprehensive
    --merge_non_CpG
    --bedGraph
    --buffer_size 10G
    --cytosine_report
    --genome_folder "$genome_folder"
    --output "$out_dir"
  )
  case "$layout" in
    paired) args+=(--paired-end) ;;
    single) args+=(--single-end) ;;
  esac

  # shellcheck disable=SC2086  # extra is intentionally word-split
  module_run "bismark_methylation_extractor" \
    bismark_methylation_extractor "${args[@]}" $extra "$bam"
}

module_bismark2report() {
  local dir="" align_rpt="" split_rpt="" mbias_rpt="" dedup_rpt="" nuc_rpt=""
  while (( $# > 0 )); do
    case "$1" in
      --dir)               dir="$2"; shift 2 ;;
      --alignment-report)  align_rpt="$2"; shift 2 ;;
      --splitting-report)  split_rpt="$2"; shift 2 ;;
      --mbias-report)      mbias_rpt="$2"; shift 2 ;;
      --dedup-report)      dedup_rpt="$2"; shift 2 ;;
      --nucleotide-report) nuc_rpt="$2"; shift 2 ;;
      *) echo "module_bismark2report: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${dir:?--dir required}"; : "${align_rpt:?--alignment-report required}"

  local args=(--dir "$dir" --alignment_report "$align_rpt")
  [[ -n "$split_rpt" ]] && args+=(--splitting_report "$split_rpt")
  [[ -n "$mbias_rpt" ]] && args+=(--mbias_report "$mbias_rpt")
  [[ -n "$dedup_rpt" ]] && args+=(--dedup_report "$dedup_rpt")
  [[ -n "$nuc_rpt"   ]] && args+=(--nucleotide_report "$nuc_rpt")

  module_run "bismark2report" bismark2report "${args[@]}"
}
