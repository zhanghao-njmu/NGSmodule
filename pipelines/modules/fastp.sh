#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/fastp.sh — adapter trimming + QC.
#
# Functions:
#   module_fastp_pe --in1 <R1> --in2 <R2> --out1 <O1> --out2 <O2>
#                   --html <H> --json <J> [--threads N] [--quality Q]
#                   [--length L] [--trim-front1 N] [--trim-front2 N]
#                   [--detect-adapter]
#   module_fastp_se --in <R> --out <O> --html <H> --json <J>
#                   [--threads N] [--quality Q] [--length L] [--trim-front N]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_fastp_pe() {
  local in1="" in2="" out1="" out2="" html="" json=""
  local threads="${threads:-4}" quality=20 length=50
  local front1=0 front2=0 detect=""
  while (( $# > 0 )); do
    case "$1" in
      --in1)         in1="$2"; shift 2 ;;
      --in2)         in2="$2"; shift 2 ;;
      --out1)        out1="$2"; shift 2 ;;
      --out2)        out2="$2"; shift 2 ;;
      --html)        html="$2"; shift 2 ;;
      --json)        json="$2"; shift 2 ;;
      --threads)     threads="$2"; shift 2 ;;
      --quality)     quality="$2"; shift 2 ;;
      --length)      length="$2"; shift 2 ;;
      --trim-front1) front1="$2"; shift 2 ;;
      --trim-front2) front2="$2"; shift 2 ;;
      --detect-adapter) detect="--detect_adapter_for_pe"; shift ;;
      *) echo "module_fastp_pe: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in1:?--in1 required}"; : "${in2:?--in2 required}"
  : "${out1:?--out1 required}"; : "${out2:?--out2 required}"

  module_run "fastp PE" \
    fastp \
      -i "$in1" -I "$in2" \
      -o "$out1" -O "$out2" \
      ${html:+--html "$html"} ${json:+--json "$json"} \
      --thread "$threads" \
      --qualified_quality_phred "$quality" \
      --length_required "$length" \
      --trim_front1 "$front1" --trim_front2 "$front2" \
      $detect
}

module_fastp_se() {
  local in="" out="" html="" json=""
  local threads="${threads:-4}" quality=20 length=50 front=0
  while (( $# > 0 )); do
    case "$1" in
      --in)         in="$2"; shift 2 ;;
      --out)        out="$2"; shift 2 ;;
      --html)       html="$2"; shift 2 ;;
      --json)       json="$2"; shift 2 ;;
      --threads)    threads="$2"; shift 2 ;;
      --quality)    quality="$2"; shift 2 ;;
      --length)     length="$2"; shift 2 ;;
      --trim-front) front="$2"; shift 2 ;;
      *) echo "module_fastp_se: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"

  module_run "fastp SE" \
    fastp \
      -i "$in" -o "$out" \
      ${html:+--html "$html"} ${json:+--json "$json"} \
      --thread "$threads" \
      --qualified_quality_phred "$quality" \
      --length_required "$length" \
      --trim_front1 "$front"
}
