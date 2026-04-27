#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/circexplorer2.sh — CIRCexplorer2 backsplice junction
# parsing + annotation. Modern alternative to find_circ / KNIFE / CIRI; uses
# STAR's chimeric output as input (no tophat2/bowtie1 dependency).
#
# Functions:
#   module_circexplorer2_parse     --aligner STAR|TopHat-Fusion|MapSplice|...
#                                  --in <chimeric_input> --out <bsj.bed>
#                                  [--paired]
#
#   module_circexplorer2_annotate  --bsj <bsj.bed> --ref <annotation.txt>
#                                  --genome <fasta> --out <circular.txt>
#
# Reference annotation format: see CIRCexplorer2 docs; built once per
# genome with `fetch_ucsc.py` or downloaded from CIRCexplorer2 site.
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_circexplorer2_parse() {
  local in="" out="" aligner="STAR"
  local paired=0
  while (( $# > 0 )); do
    case "$1" in
      --in)       in="$2"; shift 2 ;;
      --out)      out="$2"; shift 2 ;;
      --aligner)  aligner="$2"; shift 2 ;;
      --paired)   paired=1; shift ;;
      *) echo "module_circexplorer2_parse: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"

  local args=(parse -t "$aligner" -b "$out")
  (( paired )) && args+=(--pe)

  NGS_MODULE_CURRENT_TOOL=CIRCexplorer2 \
    module_run "CIRCexplorer2 parse ($aligner)" \
      CIRCexplorer2 "${args[@]}" "$in"
}

module_circexplorer2_annotate() {
  local bsj="" ref="" genome="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --bsj)     bsj="$2"; shift 2 ;;
      --ref)     ref="$2"; shift 2 ;;
      --genome)  genome="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      *) echo "module_circexplorer2_annotate: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bsj:?--bsj required}"; : "${ref:?--ref required}"
  : "${genome:?--genome required}"; : "${out:?--out required}"

  NGS_MODULE_CURRENT_TOOL=CIRCexplorer2 \
    module_run "CIRCexplorer2 annotate" \
      CIRCexplorer2 annotate -r "$ref" -g "$genome" -b "$bsj" -o "$out"
}
