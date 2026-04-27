#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/GATK_somatic_short/pipeline.sh
#
# GATK4 Mutect2 somatic short-variant pipeline. Modernises the legacy
# Analysis/GATK/GATK-somatic-short-variant.sh which used GATK3 -T MuTect2.
#
# Operating modes:
#   - Tumor-only: default (matches legacy behaviour)
#   - Tumor-normal: set $normal_sample=<other_sample_id> in ConfigFile to
#     pair the current sample's BAM (treated as tumor) with another
#     sample's BAM (treated as normal).
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

_split_known() {
  local raw="$1"
  [[ -z "$raw" ]] && return 0
  local v
  for v in $raw; do
    printf '%s\n' "$v"
  done
}

run_GATK_somatic_short_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/GATK_somatic_short"
  module_reset

  local tumor_bam
  tumor_bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ ! -f "$tumor_bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find tumor BAM: $tumor_bam (run Alignment first)"
    return 1
  fi

  # Optional normal pairing.
  local normal_bam=""
  if [[ -n "${normal_sample:-}" ]]; then
    local normal_align_dir="$work_dir/${normal_sample}/Alignment-${aligner}"
    normal_bam="$(_resolve_input_bam "${normal_sample}" "$aligner" "$normal_align_dir")"
    if [[ ! -f "$normal_bam" ]] && ! is_dry_run; then
      log_error "[$sample] paired normal BAM not found: $normal_bam"
      return 1
    fi
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi
  local genome_resolved="${genome:-${IGenome_Genome:-${IGenome_Fasta:-}}}"
  if [[ -z "$genome_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] genome (FASTA) not set and not resolvable from iGenomes"
    return 1
  fi
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"

  local prefix
  prefix="$(basename "${tumor_bam%.bam}").gatk"

  # ---- BQSR (optional) -------------------------------------------------
  local input_for_calling="$tumor_bam"
  local skip_bqsr="${skip_bqsr:-false}"

  local -a known_all=()
  while IFS= read -r v; do [[ -n "$v" ]] && known_all+=("$v"); done < <(_split_known "${known_indels[*]:-}")
  while IFS= read -r v; do [[ -n "$v" ]] && known_all+=("$v"); done < <(_split_known "${known_snps[*]:-}")

  if [[ "$skip_bqsr" != "true" ]] && (( ${#known_all[@]} > 0 )); then
    local bqsr_dir="$out_dir/BQSR"
    mkdir -p "$bqsr_dir"
    local bqsr_table="$bqsr_dir/${prefix}.BQSR.table"
    local bqsr_bam="$bqsr_dir/${prefix}.BQSR.bam"

    emit_status "[$sample] GATK somatic: BQSR BaseRecalibrator"
    local -a bqsr_args=(--bam "$tumor_bam" --ref "$genome_resolved" --out "$bqsr_table")
    local k
    for k in "${known_all[@]}"; do bqsr_args+=(--known-sites "$k"); done
    module_gatk_base_recalibrator "${bqsr_args[@]}"

    emit_status "[$sample] GATK somatic: BQSR ApplyBQSR"
    module_gatk_apply_bqsr \
      --bam "$tumor_bam" --ref "$genome_resolved" --table "$bqsr_table" --out "$bqsr_bam"
    input_for_calling="$bqsr_bam"
  else
    log_info "[$sample] BQSR skipped (skip_bqsr=$skip_bqsr, known_sites=${#known_all[@]})"
  fi

  # ---- Mutect2 ---------------------------------------------------------
  local mutect_dir="$out_dir/Mutect2"
  mkdir -p "$mutect_dir"
  local raw_vcf="$mutect_dir/${prefix}.Mutect2.vcf.gz"
  local mutect_bam="$mutect_dir/${prefix}.Mutect2.bam"

  emit_status "[$sample] GATK somatic: Mutect2"
  local -a m2_args=(
    --tumor-bam "$input_for_calling"
    --tumor-name "$sample"
    --ref "$genome_resolved"
    --out "$raw_vcf"
    --bam-out "$mutect_bam"
  )
  if [[ -n "$normal_bam" ]]; then
    m2_args+=(--normal-bam "$normal_bam" --normal-name "${normal_sample}")
  fi
  [[ -n "${germline_resource:-}"  ]] && m2_args+=(--germline-resource "$germline_resource")
  [[ -n "${panel_of_normals:-}"   ]] && m2_args+=(--panel-of-normals "$panel_of_normals")

  module_gatk_mutect2 "${m2_args[@]}"

  # ---- FilterMutectCalls ----------------------------------------------
  local filter_dir="$out_dir/FilterMutectCalls"
  mkdir -p "$filter_dir"
  local filter_vcf="$filter_dir/${prefix}.Mutect2.filter.vcf.gz"

  emit_status "[$sample] GATK somatic: FilterMutectCalls"
  module_gatk_filter_mutect_calls \
    --vcf "$raw_vcf" --ref "$genome_resolved" --out "$filter_vcf"

  if ! is_dry_run; then
    [[ -f "$raw_vcf"    ]] && emit_artifact somatic_vcf_raw      "$raw_vcf"
    [[ -f "$filter_vcf" ]] && emit_artifact somatic_vcf_filtered "$filter_vcf"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  if [[ -n "$normal_bam" ]]; then
    inputs_json="$(ngs_state_files_json --kind=aligned_bam "$tumor_bam" "$normal_bam")"
  else
    inputs_json="$(ngs_state_files_json --kind=aligned_bam "$tumor_bam")"
  fi
  outputs_json="$(ngs_state_files_json --kind=somatic_vcf_filtered "$filter_vcf")"
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","skip_bqsr":%s,"mode":"%s","normal_sample":"%s"}' \
    "$aligner" "${Deduplication:-TRUE}" "${skip_bqsr:-false}" \
    "$([[ -n "$normal_bam" ]] && echo "tumor-normal" || echo "tumor-only")" \
    "${normal_sample:-}")"

  ngs_state_stage_end "$sample" GATK_somatic_short completed \
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
  pipeline_run GATK_somatic_short
fi
