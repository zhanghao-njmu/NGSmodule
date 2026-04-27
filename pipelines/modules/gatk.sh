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

# ---------------------------------------------------------------------------
# GATK4 CNV calling toolkit. The standard germline / tumor-only flow is:
#
#   CollectReadCounts        --bam <bam> --intervals <interval_list> --out <hdf5>
#   DenoiseReadCounts        --counts <hdf5> --pon <hdf5> --std-out --denoised-out
#   ModelSegments            --denoised <tsv> --out-prefix <p> --out-dir <d>
#   CallCopyRatioSegments    --segments <cr.seg> --out <called.seg>
#   PlotModeledSegments      --segments <cr.seg> --denoised <tsv> ...
#
# Each function is a thin wrapper — the pipeline composes them.
# ---------------------------------------------------------------------------

module_gatk_collect_read_counts() {
  local bam="" intervals="" out="" extra=""
  while (( $# > 0 )); do
    case "$1" in
      --bam)        bam="$2"; shift 2 ;;
      --intervals)  intervals="$2"; shift 2 ;;
      --out)        out="$2"; shift 2 ;;
      --extra)      extra="$2"; shift 2 ;;
      *) echo "module_gatk_collect_read_counts: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${bam:?--bam required}"; : "${intervals:?--intervals required}"; : "${out:?--out required}"

  # shellcheck disable=SC2086  # extra is intentionally word-split
  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk CollectReadCounts" "$(_gatk_bin)" \
      CollectReadCounts -I "$bam" -L "$intervals" \
        --interval-merging-rule OVERLAPPING_ONLY -O "$out" $extra
}

module_gatk_denoise_read_counts() {
  local counts="" pon="" standardized_out="" denoised_out=""
  while (( $# > 0 )); do
    case "$1" in
      --counts)            counts="$2"; shift 2 ;;
      --pon)               pon="$2"; shift 2 ;;
      --standardized-out)  standardized_out="$2"; shift 2 ;;
      --denoised-out)      denoised_out="$2"; shift 2 ;;
      *) echo "module_gatk_denoise_read_counts: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${counts:?--counts required}"; : "${pon:?--pon required}"
  : "${standardized_out:?--standardized-out required}"
  : "${denoised_out:?--denoised-out required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk DenoiseReadCounts" "$(_gatk_bin)" \
      DenoiseReadCounts -I "$counts" --count-panel-of-normals "$pon" \
        --standardized-copy-ratios "$standardized_out" \
        --denoised-copy-ratios "$denoised_out"
}

module_gatk_model_segments() {
  local denoised="" out_dir="" out_prefix=""
  while (( $# > 0 )); do
    case "$1" in
      --denoised)    denoised="$2"; shift 2 ;;
      --out-dir)     out_dir="$2"; shift 2 ;;
      --out-prefix)  out_prefix="$2"; shift 2 ;;
      *) echo "module_gatk_model_segments: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${denoised:?--denoised required}"; : "${out_dir:?--out-dir required}"
  : "${out_prefix:?--out-prefix required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk ModelSegments" "$(_gatk_bin)" \
      ModelSegments \
        --denoised-copy-ratios "$denoised" \
        --output "$out_dir" --output-prefix "$out_prefix"
}

module_gatk_call_copy_ratio_segments() {
  local in="" out=""
  while (( $# > 0 )); do
    case "$1" in
      --in)  in="$2"; shift 2 ;;
      --out) out="$2"; shift 2 ;;
      *) echo "module_gatk_call_copy_ratio_segments: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk CallCopyRatioSegments" "$(_gatk_bin)" \
      CallCopyRatioSegments -I "$in" -O "$out"
}

# Optional plot generation — needs R + ggplot2 from the GATK env.
module_gatk_plot_modeled_segments() {
  local denoised="" segments="" sequence_dict="" out_dir="" out_prefix=""
  while (( $# > 0 )); do
    case "$1" in
      --denoised)        denoised="$2"; shift 2 ;;
      --segments)        segments="$2"; shift 2 ;;
      --sequence-dict)   sequence_dict="$2"; shift 2 ;;
      --out-dir)         out_dir="$2"; shift 2 ;;
      --out-prefix)      out_prefix="$2"; shift 2 ;;
      *) echo "module_gatk_plot_modeled_segments: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${denoised:?--denoised required}"; : "${segments:?--segments required}"
  : "${sequence_dict:?--sequence-dict required}"
  : "${out_dir:?--out-dir required}"; : "${out_prefix:?--out-prefix required}"

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk PlotModeledSegments" "$(_gatk_bin)" \
      PlotModeledSegments \
        --denoised-copy-ratios "$denoised" \
        --segments "$segments" \
        --sequence-dictionary "$sequence_dict" \
        --output "$out_dir" --output-prefix "$out_prefix"
}

# Project-scope helper: build a Panel of Normals from a list of read-count
# HDF5s. Caller passes one --counts <hdf5> per normal sample.
module_gatk_create_read_count_pon() {
  local out=""
  local -a counts=()
  while (( $# > 0 )); do
    case "$1" in
      --out)    out="$2"; shift 2 ;;
      --counts) counts+=("$2"); shift 2 ;;
      *) echo "module_gatk_create_read_count_pon: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${out:?--out required}"
  if (( ${#counts[@]} == 0 )); then
    echo "module_gatk_create_read_count_pon: at least one --counts required" >&2
    return 2
  fi

  local args=(CreateReadCountPanelOfNormals -O "$out")
  local c
  for c in "${counts[@]}"; do args+=(-I "$c"); done

  NGS_MODULE_CURRENT_TOOL=gatk \
    module_run "gatk CreateReadCountPanelOfNormals" "$(_gatk_bin)" "${args[@]}"
}
