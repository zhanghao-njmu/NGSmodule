#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/tabix.sh — Generic indexing for compressed VCF/BED/GFF.
#
# Functions:
#   module_tabix_index --in <bgz_file> [--preset vcf|bed|gff]
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_tabix_index() {
  local in="" preset="vcf"
  while (( $# > 0 )); do
    case "$1" in
      --in)     in="$2"; shift 2 ;;
      --preset) preset="$2"; shift 2 ;;
      *) echo "module_tabix_index: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"

  module_run "tabix -p $preset" tabix -f -p "$preset" "$in"
}
