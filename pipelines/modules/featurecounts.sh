#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/featurecounts.sh — Read summarisation per gene/exon.
#
# Wraps the `featureCounts` CLI from the subread package. The output is a
# tab-delimited count matrix (one column per BAM) plus a .summary file.
#
# Functions:
#   module_featurecounts_count --gtf <gtf> --out <counts.tab> --bam <bam>...
#                              [--threads N] [--paired] [--strand 0|1|2]
#                              [--feature-type exon] [--meta-feature gene_id]
#                              [--extra "..."]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_featurecounts_count() {
  local gtf="" out="" extra=""
  local threads="${threads:-4}"
  local paired=0
  local strand="0"
  local feature_type="exon"
  local meta_feature="gene_id"
  local -a bams=()
  while (( $# > 0 )); do
    case "$1" in
      --gtf)          gtf="$2"; shift 2 ;;
      --out)          out="$2"; shift 2 ;;
      --bam)          bams+=("$2"); shift 2 ;;
      --threads)      threads="$2"; shift 2 ;;
      --paired)       paired=1; shift ;;
      --strand)       strand="$2"; shift 2 ;;
      --feature-type) feature_type="$2"; shift 2 ;;
      --meta-feature) meta_feature="$2"; shift 2 ;;
      --extra)        extra="$2"; shift 2 ;;
      *) echo "module_featurecounts_count: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${gtf:?--gtf required}"; : "${out:?--out required}"
  if (( ${#bams[@]} == 0 )); then
    echo "module_featurecounts_count: at least one --bam required" >&2
    return 2
  fi

  local args=(
    -T "$threads"
    -a "$gtf"
    -F GTF
    -t "$feature_type"
    -g "$meta_feature"
    -s "$strand"
    -o "$out"
  )
  (( paired )) && args+=(-p --countReadPairs)

  # shellcheck disable=SC2086  # extra is intentionally word-split
  module_run "featureCounts" \
    featureCounts "${args[@]}" $extra "${bams[@]}"
}
