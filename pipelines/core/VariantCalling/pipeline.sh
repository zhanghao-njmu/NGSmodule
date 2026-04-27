#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/VariantCalling/pipeline.sh
#
# bcftools mpileup → call → split (snps|indels) → MQ-filter → concat.
# Mirrors the legacy Analysis/Samtools/Samtools-variant-calling.sh logic
# but composed from modular bcftools/tabix wrappers with full provenance.
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
module_load bcftools
module_load tabix

# ---------------------------------------------------------------------------
# Resolve which BAM to call against based on Deduplication state.
# ---------------------------------------------------------------------------
_resolve_input_bam() {
  local sample="$1" aligner="$2" align_dir="$3"
  case "${Deduplication:-TRUE}" in
    FALSE) echo "$align_dir/${sample}.${aligner}.markdup.bam" ;;
    TRUE)  echo "$align_dir/${sample}.${aligner}.dedup.bam" ;;
    *)     echo "$align_dir/${sample}.${aligner}.sorted.bam" ;;
  esac
}

# ---------------------------------------------------------------------------
# Per-sample worker.
# ---------------------------------------------------------------------------
run_VariantCalling_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/VariantCalling"
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
  if [[ -n "$genome_resolved" ]] && [[ ! -f "$genome_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] genome FASTA not found: $genome_resolved"
    return 1
  fi

  local prefix
  prefix="$(basename "${bam%.bam}").bcftools"
  local raw_vcf="$out_dir/Samtools/${prefix}.vcf.gz"
  local raw_stats="${raw_vcf}.stats"
  mkdir -p "$out_dir/Samtools" "$out_dir/VariantFiltration"

  # ---- Mpileup + call ---------------------------------------------------
  emit_status "[$sample] VariantCalling: bcftools mpileup | call"
  module_bcftools_mpileup_call \
    --bam "$bam" --ref "${genome_resolved:-<genome>}" \
    --vcf "$raw_vcf" --threads "${threads:-4}"

  module_tabix_index --in "$raw_vcf" --preset vcf

  emit_status "[$sample] VariantCalling: bcftools stats"
  module_bcftools_stats \
    --vcf "$raw_vcf" --ref "${genome_resolved:-<genome>}" --out "$raw_stats"

  # ---- Split → filter → concat ----------------------------------------
  local snps_raw="$out_dir/VariantFiltration/${prefix}.snps.vcf.gz"
  local indels_raw="$out_dir/VariantFiltration/${prefix}.indels.vcf.gz"
  local snps_filt="$out_dir/VariantFiltration/${prefix}.snps.filter.vcf.gz"
  local indels_filt="$out_dir/VariantFiltration/${prefix}.indels.filter.vcf.gz"
  local final_vcf="$out_dir/VariantFiltration/${prefix}.filter.vcf.gz"

  emit_status "[$sample] VariantCalling: split snps + indels"
  module_bcftools_view_filter --in "$raw_vcf" --out "$snps_raw"   --type snps   --threads "${threads:-4}"
  module_bcftools_view_filter --in "$raw_vcf" --out "$indels_raw" --type indels --threads "${threads:-4}"

  emit_status "[$sample] VariantCalling: MQ filter"
  module_bcftools_filter --in "$snps_raw"   --out "$snps_filt"   --expr "MQ < ${snp_mq:-40}"   --threads "${threads:-4}"
  module_bcftools_filter --in "$indels_raw" --out "$indels_filt" --expr "MQ < ${indel_mq:-20}" --threads "${threads:-4}"

  module_bcftools_index --vcf "$snps_filt"
  module_bcftools_index --vcf "$indels_filt"

  emit_status "[$sample] VariantCalling: concat filtered"
  module_bcftools_concat --out "$final_vcf" \
    --in "$snps_filt" --in "$indels_filt" --threads "${threads:-4}"

  if ! is_dry_run; then
    emit_artifact variant_vcf_raw      "$raw_vcf"
    emit_artifact variant_vcf_filtered "$final_vcf"
    [[ -f "$raw_stats" ]] && emit_artifact variant_stats "$raw_stats"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  outputs_json="$(ngs_state_files_json --kind=variant_vcf_filtered "$final_vcf")"
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","snp_mq":%s,"indel_mq":%s}' \
    "$aligner" "${Deduplication:-TRUE}" "${snp_mq:-40}" "${indel_mq:-20}")"

  ngs_state_stage_end "$sample" VariantCalling completed \
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
  pipeline_run VariantCalling
fi
