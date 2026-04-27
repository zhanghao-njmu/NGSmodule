#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/GATK_germline_short/pipeline.sh
#
# GATK4 germline short-variant pipeline. Modernises the legacy
# Analysis/GATK/GATK-germline-short-variant.sh which used GATK3 syntax.
#
# Stages per sample:
#   1. BQSR    — BaseRecalibrator + ApplyBQSR (skipped if no known sites)
#   2. Calling — HaplotypeCaller (GVCF) + GenotypeGVCFs
#   3. Filter  — split SNP/INDEL → VariantFiltration → MergeVcfs
###############################################################################

set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

# shellcheck source=../../lib/orchestrator.sh
source "$_NGS_LIB_DIR/orchestrator.sh"

_NGS_MODULES_DIR="$(cd "$_PIPELINE_DIR/../../modules" && pwd)"
# shellcheck source=../../modules/_lib.sh
source "$_NGS_MODULES_DIR/_lib.sh"
module_load gatk

_resolve_input_bam() {
  local sample="$1" aligner="$2" align_dir="$3"
  case "${Deduplication:-TRUE}" in
    FALSE) echo "$align_dir/${sample}.${aligner}.markdup.bam" ;;
    TRUE)  echo "$align_dir/${sample}.${aligner}.dedup.bam" ;;
    *)     echo "$align_dir/${sample}.${aligner}.sorted.bam" ;;
  esac
}

# ---------------------------------------------------------------------------
# Read whitespace-separated VCF paths from $known_indels / $known_snps
# (config-supplied) into a flat list. Empty list returned if unset.
# ---------------------------------------------------------------------------
_split_known() {
  local raw="$1"
  [[ -z "$raw" ]] && return 0
  local v
  for v in $raw; do
    printf '%s\n' "$v"
  done
}

run_GATK_germline_short_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/GATK_germline_short"
  module_reset

  local bam
  bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ ! -f "$bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find input BAM: $bam (run Alignment first)"
    return 1
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi

  local genome_resolved="${genome:-${IGenome_Genome:-${IGenome_Fasta:-}}}"
  if [[ -z "$genome_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] genome (FASTA) not set and not resolvable from iGenomes"
    return 1
  fi
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"   # dry-run placeholder

  local prefix
  prefix="$(basename "${bam%.bam}").gatk"

  # ---- Step 1: BQSR ----------------------------------------------------
  local input_for_calling="$bam"
  local skip_bqsr="${skip_bqsr:-false}"

  # Aggregate known sites from the optional config arrays/strings.
  local -a known_all=()
  while IFS= read -r v; do [[ -n "$v" ]] && known_all+=("$v"); done < <(_split_known "${known_indels[*]:-}")
  while IFS= read -r v; do [[ -n "$v" ]] && known_all+=("$v"); done < <(_split_known "${known_snps[*]:-}")

  if [[ "$skip_bqsr" != "true" ]] && (( ${#known_all[@]} > 0 )); then
    local bqsr_dir="$out_dir/BQSR"
    mkdir -p "$bqsr_dir"
    local bqsr_table="$bqsr_dir/${prefix}.BQSR.table"
    local bqsr_bam="$bqsr_dir/${prefix}.BQSR.bam"

    emit_status "[$sample] GATK germline: BQSR BaseRecalibrator"
    local -a bqsr_args=(--bam "$bam" --ref "$genome_resolved" --out "$bqsr_table")
    local k
    for k in "${known_all[@]}"; do bqsr_args+=(--known-sites "$k"); done
    module_gatk_base_recalibrator "${bqsr_args[@]}"

    emit_status "[$sample] GATK germline: BQSR ApplyBQSR"
    module_gatk_apply_bqsr \
      --bam "$bam" --ref "$genome_resolved" --table "$bqsr_table" --out "$bqsr_bam"
    input_for_calling="$bqsr_bam"
  else
    log_info "[$sample] BQSR skipped (skip_bqsr=$skip_bqsr, known_sites=${#known_all[@]})"
  fi

  # ---- Step 2: HaplotypeCaller + GenotypeGVCFs ------------------------
  local hc_dir="$out_dir/HaplotypeCaller"
  mkdir -p "$hc_dir"
  local gvcf="$hc_dir/${prefix}.HaplotypeCaller.gvcf.gz"
  local raw_vcf="$hc_dir/${prefix}.HaplotypeCaller.vcf.gz"
  local hc_bam="$hc_dir/${prefix}.HaplotypeCaller.bam"

  emit_status "[$sample] GATK germline: HaplotypeCaller (GVCF)"
  module_gatk_haplotype_caller \
    --bam "$input_for_calling" --ref "$genome_resolved" \
    --out "$gvcf" --bam-out "$hc_bam" --gvcf

  emit_status "[$sample] GATK germline: GenotypeGVCFs"
  module_gatk_genotype_gvcfs \
    --vcf "$gvcf" --ref "$genome_resolved" --out "$raw_vcf"

  # ---- Step 3: Hard-filter SNPs + INDELs separately, merge ------------
  local vf_dir="$out_dir/VariantFiltration"
  mkdir -p "$vf_dir"
  local snps_raw="$vf_dir/${prefix}.snps.vcf.gz"
  local indels_raw="$vf_dir/${prefix}.indels.vcf.gz"
  local snps_filt="$vf_dir/${prefix}.snps.filter.vcf.gz"
  local indels_filt="$vf_dir/${prefix}.indels.filter.vcf.gz"
  local final_vcf="$vf_dir/${prefix}.HaplotypeCaller.filter.vcf.gz"

  emit_status "[$sample] GATK germline: split SNP/INDEL"
  module_gatk_select_variants \
    --vcf "$raw_vcf" --ref "$genome_resolved" --out "$snps_raw"   --select-type SNP
  module_gatk_select_variants \
    --vcf "$raw_vcf" --ref "$genome_resolved" --out "$indels_raw" --select-type INDEL

  emit_status "[$sample] GATK germline: VariantFiltration"
  module_gatk_variant_filtration \
    --vcf "$snps_raw" --ref "$genome_resolved" --out "$snps_filt" \
    --filter-expression "${snp_filter:?snp_filter not set}" \
    --filter-name "snp_hard_filter"
  module_gatk_variant_filtration \
    --vcf "$indels_raw" --ref "$genome_resolved" --out "$indels_filt" \
    --filter-expression "${indel_filter:?indel_filter not set}" \
    --filter-name "indel_hard_filter"

  emit_status "[$sample] GATK germline: MergeVcfs"
  module_gatk_merge_vcfs --out "$final_vcf" --in "$snps_filt" --in "$indels_filt"

  if ! is_dry_run; then
    emit_artifact gvcf                 "$gvcf"
    emit_artifact variant_vcf          "$raw_vcf"
    emit_artifact variant_vcf_filtered "$final_vcf"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  outputs_json="$(ngs_state_files_json --kind=variant_vcf_filtered "$final_vcf")"
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","skip_bqsr":%s,"known_sites_count":%d}' \
    "$aligner" "${Deduplication:-TRUE}" "${skip_bqsr:-false}" "${#known_all[@]}")"

  ngs_state_stage_end "$sample" GATK_germline_short completed \
    --params  "$params_json" \
    --tools   "$tools_json" \
    --inputs  "$inputs_json" \
    --outputs "$outputs_json"

  return 0
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  if [[ -n "${ConfigFile:-}" ]] && [[ -f "${ConfigFile}" ]]; then
    # shellcheck source=/dev/null
    source "$ConfigFile"
  fi
  if [[ -n "${SampleInfoFile:-}" ]] && [[ -f "${SampleInfoFile}" ]]; then
    ngs_load_sample_info "$SampleInfoFile"
  fi
  if [[ -n "${rawdata_dir:-}" ]] && [[ -n "${work_dir:-}" ]]; then
    if [[ ! -d "$work_dir" ]] || [[ -z "$(ls -A "$work_dir" 2>/dev/null)" ]]; then
      ngs_build_workdir "$rawdata_dir" "$work_dir"
    fi
  fi
  pipeline_run GATK_germline_short
fi
