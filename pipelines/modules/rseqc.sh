#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/rseqc.sh — Quality assessment of aligned RNA-seq.
#
# RSeQC is a collection of Python scripts; we expose the ones the legacy
# postAlignmentQC pipeline used. Each function writes to ${out-prefix}.* in
# the current working directory or an explicit --out-prefix.
#
# Functions:
#   module_rseqc_bam_stat            --bam <bam> --out <txt>
#   module_rseqc_infer_experiment    --bam <bam> --bed <bed> --out <txt>
#   module_rseqc_inner_distance      --bam <bam> --bed <bed> --out-prefix <p>
#   module_rseqc_read_distribution   --bam <bam> --bed <bed> --out <txt>
#   module_rseqc_read_duplication    --bam <bam> --out-prefix <p>
#   module_rseqc_read_gc             --bam <bam> --out-prefix <p>
#   module_rseqc_genebody_coverage   --bam <bam> --bed <bed> --out-prefix <p>
#   module_rseqc_junction_annotation --bam <bam> --bed <bed> --out-prefix <p>
#   module_rseqc_junction_saturation --bam <bam> --bed <bed> --out-prefix <p>
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

_rseqc_parse() {
  local _name="$1"; shift
  bam=""; bed=""; out=""; out_prefix=""
  while (( $# > 0 )); do
    case "$1" in
      --bam)        bam="$2"; shift 2 ;;
      --bed)        bed="$2"; shift 2 ;;
      --out)        out="$2"; shift 2 ;;
      --out-prefix) out_prefix="$2"; shift 2 ;;
      *) echo "$_name: unknown arg $1" >&2; return 2 ;;
    esac
  done
}

module_rseqc_bam_stat() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_bam_stat "$@" || return $?
  : "${bam:?--bam required}"; : "${out:?--out required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc bam_stat" bash -c "bam_stat.py -i '$bam' > '$out' 2>&1"
}

module_rseqc_infer_experiment() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_infer_experiment "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out:?--out required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc infer_experiment" bash -c \
      "infer_experiment.py -r '$bed' -i '$bam' > '$out' 2>&1"
}

module_rseqc_inner_distance() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_inner_distance "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc inner_distance" \
      inner_distance.py -r "$bed" -i "$bam" -o "$out_prefix"
}

module_rseqc_read_distribution() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_read_distribution "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out:?--out required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc read_distribution" bash -c \
      "read_distribution.py -r '$bed' -i '$bam' > '$out' 2>&1"
}

module_rseqc_read_duplication() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_read_duplication "$@" || return $?
  : "${bam:?--bam required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc read_duplication" \
      read_duplication.py -i "$bam" -o "$out_prefix"
}

module_rseqc_read_gc() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_read_gc "$@" || return $?
  : "${bam:?--bam required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc read_GC" \
      read_GC.py -i "$bam" -o "$out_prefix"
}

module_rseqc_genebody_coverage() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_genebody_coverage "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc geneBody_coverage" \
      geneBody_coverage.py -r "$bed" -i "$bam" -o "$out_prefix"
}

module_rseqc_junction_annotation() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_junction_annotation "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc junction_annotation" \
      junction_annotation.py -r "$bed" -i "$bam" -o "$out_prefix"
}

module_rseqc_junction_saturation() {
  local bam bed out out_prefix
  _rseqc_parse module_rseqc_junction_saturation "$@" || return $?
  : "${bam:?--bam required}"; : "${bed:?--bed required}"; : "${out_prefix:?--out-prefix required}"
  NGS_MODULE_CURRENT_TOOL=rseqc \
    module_run "rseqc junction_saturation" \
      junction_saturation.py -r "$bed" -i "$bam" -o "$out_prefix"
}
