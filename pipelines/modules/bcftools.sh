#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/bcftools.sh — VCF/BCF manipulation.
#
# Functions (each maps 1:1 to a bcftools subcommand):
#   module_bcftools_mpileup_call --bam <bam> --ref <fa> --vcf <out.vcf.gz>
#                                [--threads N] [--max-depth D]
#   module_bcftools_view_filter  --in <vcf.gz> --out <vcf.gz> --type snps|indels
#                                [--threads N]
#   module_bcftools_filter       --in <vcf.gz> --out <vcf.gz> --expr <expr>
#                                [--threads N]
#   module_bcftools_concat       --out <vcf.gz> --in <vcf.gz>... [--threads N]
#   module_bcftools_index        --vcf <vcf.gz>
#   module_bcftools_stats        --vcf <vcf.gz> --ref <fa> --out <txt>
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

# Pipe `bcftools mpileup | bcftools call -vmO z` into one bgzip-compressed VCF.
module_bcftools_mpileup_call() {
  local bam="" ref="" vcf=""
  local threads="${threads:-4}"
  local max_depth=100000000
  while (( $# > 0 )); do
    case "$1" in
      --bam)       bam="$2"; shift 2 ;;
      --ref)       ref="$2"; shift 2 ;;
      --vcf)       vcf="$2"; shift 2 ;;
      --threads)   threads="$2"; shift 2 ;;
      --max-depth) max_depth="$2"; shift 2 ;;
      *) echo "module_bcftools_mpileup_call: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${ref:?--ref required}"; : "${vcf:?--vcf required}"

  NGS_MODULE_CURRENT_TOOL=bcftools \
    module_run "bcftools mpileup | call" bash -c \
      "bcftools mpileup --threads '$threads' -d '$max_depth' -Ou -f '$ref' '$bam' \
       | bcftools call --threads '$threads' -vmO z -o '$vcf'"
}

module_bcftools_view_filter() {
  local in="" out="" type=""
  local threads="${threads:-4}"
  while (( $# > 0 )); do
    case "$1" in
      --in)      in="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      --type)    type="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *) echo "module_bcftools_view_filter: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"; : "${type:?--type required}"

  module_run "bcftools view --types $type" \
    bcftools view --threads "$threads" --types "$type" \
                  --output-type z --output-file "$out" "$in"
}

module_bcftools_filter() {
  local in="" out="" expr=""
  local threads="${threads:-4}"
  while (( $# > 0 )); do
    case "$1" in
      --in)      in="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      --expr)    expr="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *) echo "module_bcftools_filter: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"; : "${expr:?--expr required}"

  module_run "bcftools filter -e '$expr'" \
    bcftools filter --threads "$threads" -e "$expr" \
                    --output-type z --output "$out" "$in"
}

module_bcftools_concat() {
  local out=""
  local threads="${threads:-4}"
  local -a ins=()
  while (( $# > 0 )); do
    case "$1" in
      --out)     out="$2"; shift 2 ;;
      --in)      ins+=("$2"); shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *) echo "module_bcftools_concat: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${out:?--out required}"
  if (( ${#ins[@]} == 0 )); then
    echo "module_bcftools_concat: at least one --in required" >&2
    return 2
  fi

  module_run "bcftools concat" \
    bcftools concat --threads "$threads" -a \
                    --output-type z --output "$out" "${ins[@]}"
}

module_bcftools_index() {
  local vcf=""
  while (( $# > 0 )); do
    case "$1" in
      --vcf) vcf="$2"; shift 2 ;;
      *) echo "module_bcftools_index: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"

  module_run "bcftools index" bcftools index -f "$vcf"
}

module_bcftools_stats() {
  local vcf="" ref="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --vcf) vcf="$2"; shift 2 ;;
      --ref) ref="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      *) echo "module_bcftools_stats: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"; : "${out:?--out required}"

  local args=(-s -)
  [[ -n "$ref" ]] && args+=(-F "$ref")
  NGS_MODULE_CURRENT_TOOL=bcftools \
    module_run "bcftools stats" bash -c \
      "bcftools stats ${args[*]} '$vcf' > '$out'"
}
