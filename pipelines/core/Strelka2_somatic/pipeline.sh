#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/Strelka2_somatic/pipeline.sh
#
# Strelka2 tumor-normal pair somatic variant calling. Unlike the legacy
# Analysis/Strelka2/Strelka2-somatic-short-variant.sh (which forgot to
# wire up --normalBam), this implementation requires the matched normal
# explicitly via $normal_sample. Strelka2's somatic mode is invalid
# without a normal; we surface that as a hard error.
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
module_load strelka2
module_load bcftools

_resolve_input_bam() {
  local sample="$1" aligner="$2" align_dir="$3"
  case "${Deduplication:-TRUE}" in
    FALSE) echo "$align_dir/${sample}.${aligner}.markdup.bam" ;;
    TRUE)  echo "$align_dir/${sample}.${aligner}.dedup.bam" ;;
    *)     echo "$align_dir/${sample}.${aligner}.sorted.bam" ;;
  esac
}

run_Strelka2_somatic_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/Strelka2_somatic"
  module_reset

  if [[ -z "${normal_sample:-}" ]]; then
    log_error "[$sample] normal_sample is required for somatic calling (set in ConfigFile)"
    return 1
  fi

  local tumor_bam normal_bam
  tumor_bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ ! -f "$tumor_bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find tumor BAM: $tumor_bam (run Alignment first)"
    return 1
  fi

  local normal_align_dir="$work_dir/${normal_sample}/Alignment-${aligner}"
  normal_bam="$(_resolve_input_bam "${normal_sample}" "$aligner" "$normal_align_dir")"
  if [[ ! -f "$normal_bam" ]] && ! is_dry_run; then
    log_error "[$sample] paired normal BAM not found: $normal_bam"
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
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"

  local prefix
  prefix="$(basename "${tumor_bam%.bam}").Strelka2"

  local rundir="$out_dir/Strelka2"
  if ! is_dry_run; then
    rm -rf "$rundir"
  fi
  mkdir -p "$out_dir"

  emit_status "[$sample] Strelka2 somatic: configure + run (tumor=$sample, normal=$normal_sample)"
  local -a args=(
    --normal-bam "$normal_bam"
    --tumor-bam "$tumor_bam"
    --ref "$genome_resolved"
    --rundir "$rundir"
    --threads "${threads:-4}"
  )
  [[ "${exome:-FALSE}"    == "TRUE" ]] && args+=(--exome)
  [[ "${targeted:-FALSE}" == "TRUE" ]] && args+=(--targeted)

  module_strelka2_somatic "${args[@]}"

  # Strelka2 somatic produces two VCFs: somatic.snvs.vcf.gz and somatic.indels.vcf.gz
  local final_snvs="$out_dir/${prefix}.somatic.snvs.vcf.gz"
  local final_indels="$out_dir/${prefix}.somatic.indels.vcf.gz"

  if ! is_dry_run; then
    if [[ -f "$rundir/results/variants/somatic.snvs.vcf.gz" ]]; then
      cp "$rundir/results/variants/somatic.snvs.vcf.gz"       "$final_snvs"
      cp "$rundir/results/variants/somatic.snvs.vcf.gz.tbi"   "${final_snvs}.tbi"
      cp "$rundir/results/variants/somatic.indels.vcf.gz"     "$final_indels"
      cp "$rundir/results/variants/somatic.indels.vcf.gz.tbi" "${final_indels}.tbi"
    fi
  fi

  emit_status "[$sample] Strelka2 somatic: bcftools stats"
  module_bcftools_stats --vcf "$final_snvs"   --ref "$genome_resolved" --out "${final_snvs}.stats"
  module_bcftools_stats --vcf "$final_indels" --ref "$genome_resolved" --out "${final_indels}.stats"

  if ! is_dry_run; then
    [[ -f "$final_snvs"   ]] && emit_artifact somatic_snvs_vcf   "$final_snvs"
    [[ -f "$final_indels" ]] && emit_artifact somatic_indels_vcf "$final_indels"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$tumor_bam" "$normal_bam")"
  outputs_json="$(ngs_state_files_json --kind=somatic_snvs_vcf "$final_snvs" "$final_indels")"
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","exome":"%s","targeted":"%s","normal_sample":"%s"}' \
    "$aligner" "${Deduplication:-TRUE}" "${exome:-FALSE}" "${targeted:-FALSE}" "${normal_sample}")"

  ngs_state_stage_end "$sample" Strelka2_somatic completed \
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
  pipeline_run Strelka2_somatic
fi
