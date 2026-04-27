#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/strelka2.sh — Strelka2 small-variant caller.
#
# Strelka2's CLI is two-step: first configureStrelka{Germline,Somatic}Workflow.py
# generates a runDir, then runDir/runWorkflow.py executes it. This module
# wraps both into a single function call per use case.
#
# Functions:
#   module_strelka2_germline --bam <bam> --ref <fa> --rundir <dir>
#                            [--exome] [--targeted] [--threads N]
#
#   module_strelka2_somatic  --normal-bam <bam> --tumor-bam <bam> --ref <fa>
#                            --rundir <dir> [--exome] [--targeted] [--threads N]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_strelka2_germline() {
  local bam="" ref="" rundir=""
  local threads="${threads:-4}"
  local exome=0 targeted=0
  while (( $# > 0 )); do
    case "$1" in
      --bam)      bam="$2"; shift 2 ;;
      --ref)      ref="$2"; shift 2 ;;
      --rundir)   rundir="$2"; shift 2 ;;
      --threads)  threads="$2"; shift 2 ;;
      --exome)    exome=1; shift ;;
      --targeted) targeted=1; shift ;;
      *) echo "module_strelka2_germline: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${ref:?--ref required}"; : "${rundir:?--rundir required}"

  local cfg_args=(--bam "$bam" --referenceFasta "$ref" --runDir "$rundir")
  (( exome ))    && cfg_args+=(--exome)
  (( targeted )) && cfg_args+=(--targeted)

  NGS_MODULE_CURRENT_TOOL=strelka2 \
    module_run "Strelka2 germline configure" \
      configureStrelkaGermlineWorkflow.py "${cfg_args[@]}"

  NGS_MODULE_CURRENT_TOOL=strelka2 \
    module_run "Strelka2 germline runWorkflow" \
      "$rundir/runWorkflow.py" -m local -j "$threads"
}

module_strelka2_somatic() {
  local normal_bam="" tumor_bam="" ref="" rundir=""
  local threads="${threads:-4}"
  local exome=0 targeted=0
  while (( $# > 0 )); do
    case "$1" in
      --normal-bam) normal_bam="$2"; shift 2 ;;
      --tumor-bam)  tumor_bam="$2"; shift 2 ;;
      --ref)        ref="$2"; shift 2 ;;
      --rundir)     rundir="$2"; shift 2 ;;
      --threads)    threads="$2"; shift 2 ;;
      --exome)      exome=1; shift ;;
      --targeted)   targeted=1; shift ;;
      *) echo "module_strelka2_somatic: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${normal_bam:?--normal-bam required}"; : "${tumor_bam:?--tumor-bam required}"
  : "${ref:?--ref required}"; : "${rundir:?--rundir required}"

  local cfg_args=(
    --normalBam "$normal_bam" --tumorBam "$tumor_bam"
    --referenceFasta "$ref" --runDir "$rundir"
  )
  (( exome ))    && cfg_args+=(--exome)
  (( targeted )) && cfg_args+=(--targeted)

  NGS_MODULE_CURRENT_TOOL=strelka2 \
    module_run "Strelka2 somatic configure" \
      configureStrelkaSomaticWorkflow.py "${cfg_args[@]}"

  NGS_MODULE_CURRENT_TOOL=strelka2 \
    module_run "Strelka2 somatic runWorkflow" \
      "$rundir/runWorkflow.py" -m local -j "$threads"
}
