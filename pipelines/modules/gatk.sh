#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/gatk.sh — GATK4 wrappers (Genome Analysis Toolkit).
#
# The legacy NGSmodule scripts used GATK3 with -T <walker> syntax. This
# module modernises to GATK4 (`gatk <Tool>`); the call sites in
# pipelines/core/<...> are responsible for adapting argument differences
# (e.g. -knownSites → --known-sites, -nct → --tmp-dir + parallelism flags).
#
# Functions (1:1 with the GATK4 Tool name):
#   module_gatk_base_recalibrator   --bam --ref --out [--known-sites <vcf>]...
#   module_gatk_apply_bqsr          --bam --ref --table --out
#   module_gatk_haplotype_caller    --bam --ref --out [--gvcf] [--bam-out <bam>]
#   module_gatk_genotype_gvcfs      --vcf --ref --out
#   module_gatk_select_variants     --vcf --ref --out --select-type SNP|INDEL
#   module_gatk_variant_filtration  --vcf --ref --out --filter-expression --filter-name
#   module_gatk_merge_vcfs          --out --in <vcf>...
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

# Some sites pin GATK to a Java jar. The wrapper honours $GATK_BIN if set;
# otherwise falls back to `gatk` on PATH.
_gatk_bin() { echo "${GATK_BIN:-gatk}"; }

module_gatk_base_recalibrator() {
  local bam="" ref="" out=""
  local -a known=()
  while (( $# > 0 )); do
    case "$1" in
      --bam)         bam="$2"; shift 2 ;;
      --ref)         ref="$2"; shift 2 ;;
      --out)         out="$2"; shift 2 ;;
      --known-sites) known+=("$2"); shift 2 ;;
      *) echo "module_gatk_base_recalibrator: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${ref:?--ref required}"; : "${out:?--out required}"

  local args=(BaseRecalibrator -R "$ref" -I "$bam" -O "$out")
  local k
  for k in "${known[@]}"; do args+=(--known-sites "$k"); done

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk BaseRecalibrator" "$(_gatk_bin)" "${args[@]}"
}

module_gatk_apply_bqsr() {
  local bam="" ref="" table="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --bam)   bam="$2"; shift 2 ;;
      --ref)   ref="$2"; shift 2 ;;
      --table) table="$2"; shift 2 ;;
      --out)   out="$2"; shift 2 ;;
      *) echo "module_gatk_apply_bqsr: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${ref:?--ref required}"
  : "${table:?--table required}"; : "${out:?--out required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk ApplyBQSR" "$(_gatk_bin)" \
      ApplyBQSR -R "$ref" -I "$bam" --bqsr-recal-file "$table" -O "$out"
}

module_gatk_haplotype_caller() {
  local bam="" ref="" out="" bam_out="" gvcf=0
  while (( $# > 0 )); do
    case "$1" in
      --bam)     bam="$2"; shift 2 ;;
      --ref)     ref="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      --bam-out) bam_out="$2"; shift 2 ;;
      --gvcf)    gvcf=1; shift ;;
      *) echo "module_gatk_haplotype_caller: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${ref:?--ref required}"; : "${out:?--out required}"

  local args=(HaplotypeCaller -R "$ref" -I "$bam" -O "$out")
  (( gvcf ))           && args+=(--emit-ref-confidence GVCF)
  [[ -n "$bam_out" ]]  && args+=(--bam-output "$bam_out")

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk HaplotypeCaller" "$(_gatk_bin)" "${args[@]}"
}

module_gatk_genotype_gvcfs() {
  local vcf="" ref="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --vcf) vcf="$2"; shift 2 ;;
      --ref) ref="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      *) echo "module_gatk_genotype_gvcfs: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"; : "${ref:?--ref required}"; : "${out:?--out required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk GenotypeGVCFs" "$(_gatk_bin)" \
      GenotypeGVCFs -R "$ref" -V "$vcf" -O "$out"
}

module_gatk_select_variants() {
  local vcf="" ref="" out="" select_type=""
  while (( $# > 0 )); do
    case "$1" in
      --vcf)         vcf="$2"; shift 2 ;;
      --ref)         ref="$2"; shift 2 ;;
      --out)         out="$2"; shift 2 ;;
      --select-type) select_type="$2"; shift 2 ;;
      *) echo "module_gatk_select_variants: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"; : "${ref:?--ref required}"
  : "${out:?--out required}"; : "${select_type:?--select-type required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk SelectVariants ($select_type)" "$(_gatk_bin)" \
      SelectVariants -R "$ref" -V "$vcf" --select-type-to-include "$select_type" -O "$out"
}

module_gatk_variant_filtration() {
  local vcf="" ref="" out="" expr="" name="VariantFiltration"
  while (( $# > 0 )); do
    case "$1" in
      --vcf)               vcf="$2"; shift 2 ;;
      --ref)               ref="$2"; shift 2 ;;
      --out)               out="$2"; shift 2 ;;
      --filter-expression) expr="$2"; shift 2 ;;
      --filter-name)       name="$2"; shift 2 ;;
      *) echo "module_gatk_variant_filtration: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"; : "${ref:?--ref required}"
  : "${out:?--out required}"; : "${expr:?--filter-expression required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk VariantFiltration ($name)" "$(_gatk_bin)" \
      VariantFiltration -R "$ref" -V "$vcf" \
        --filter-expression "$expr" --filter-name "$name" -O "$out"
}

module_gatk_merge_vcfs() {
  local out=""
  local -a ins=()
  while (( $# > 0 )); do
    case "$1" in
      --out) out="$2"; shift 2 ;;
      --in)  ins+=("$2"); shift 2 ;;
      *) echo "module_gatk_merge_vcfs: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${out:?--out required}"
  if (( ${#ins[@]} == 0 )); then
    echo "module_gatk_merge_vcfs: at least one --in required" >&2
    return 2
  fi

  local args=(MergeVcfs -O "$out")
  local v
  for v in "${ins[@]}"; do args+=(-I "$v"); done

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk MergeVcfs" "$(_gatk_bin)" "${args[@]}"
}

# ---------------------------------------------------------------------------
# Somatic variant calling. Tumor-only mode is the default (matches the legacy
# Analysis/GATK/GATK-somatic-short-variant.sh which used GATK3 -T MuTect2).
# Pass --normal-bam to enable proper tumor-normal pair calling.
# ---------------------------------------------------------------------------
module_gatk_mutect2() {
  local tumor_bam="" normal_bam="" tumor_name="" normal_name=""
  local ref="" out="" bam_out="" germline_resource="" panel_of_normals=""
  while (( $# > 0 )); do
    case "$1" in
      --tumor-bam)         tumor_bam="$2"; shift 2 ;;
      --normal-bam)        normal_bam="$2"; shift 2 ;;
      --tumor-name)        tumor_name="$2"; shift 2 ;;
      --normal-name)       normal_name="$2"; shift 2 ;;
      --ref)               ref="$2"; shift 2 ;;
      --out)               out="$2"; shift 2 ;;
      --bam-out)           bam_out="$2"; shift 2 ;;
      --germline-resource) germline_resource="$2"; shift 2 ;;
      --panel-of-normals)  panel_of_normals="$2"; shift 2 ;;
      *) echo "module_gatk_mutect2: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${tumor_bam:?--tumor-bam required}"; : "${ref:?--ref required}"; : "${out:?--out required}"

  local args=(Mutect2 -R "$ref" -I "$tumor_bam" -O "$out")
  [[ -n "$tumor_name"        ]] && args+=(-tumor "$tumor_name")
  [[ -n "$normal_bam"        ]] && args+=(-I "$normal_bam")
  [[ -n "$normal_name"       ]] && args+=(-normal "$normal_name")
  [[ -n "$bam_out"           ]] && args+=(--bam-output "$bam_out")
  [[ -n "$germline_resource" ]] && args+=(--germline-resource "$germline_resource")
  [[ -n "$panel_of_normals"  ]] && args+=(--panel-of-normals "$panel_of_normals")

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk Mutect2" "$(_gatk_bin)" "${args[@]}"
}

# Mutect2 emits a stats file alongside the VCF that FilterMutectCalls needs;
# we rely on GATK's default stats path (<vcf>.stats) — caller doesn't pass it.
module_gatk_filter_mutect_calls() {
  local vcf="" ref="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --vcf) vcf="$2"; shift 2 ;;
      --ref) ref="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      *) echo "module_gatk_filter_mutect_calls: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${vcf:?--vcf required}"; : "${ref:?--ref required}"; : "${out:?--out required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk FilterMutectCalls" "$(_gatk_bin)" \
      FilterMutectCalls -R "$ref" -V "$vcf" -O "$out"
}
