#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/Strelka2_germline/pipeline.sh
#
# Strelka2 germline small-variant calling. Strelka2 produces its results
# under <rundir>/results/variants/variants.vcf.gz which we copy out to a
# stable filename matching our naming conventions.
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

run_Strelka2_germline_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/Strelka2_germline"
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
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"

  local prefix
  prefix="$(basename "${bam%.bam}").Strelka2"

  local rundir="$out_dir/Strelka2"
  if ! is_dry_run; then
    rm -rf "$rundir"
  fi
  mkdir -p "$out_dir"

  emit_status "[$sample] Strelka2 germline: configure + run"
  local -a args=(
    --bam "$bam"
    --ref "$genome_resolved"
    --rundir "$rundir"
    --threads "${threads:-4}"
  )
  [[ "${exome:-FALSE}"    == "TRUE" ]] && args+=(--exome)
  [[ "${targeted:-FALSE}" == "TRUE" ]] && args+=(--targeted)

  module_strelka2_germline "${args[@]}"

  # Strelka2 writes results/variants/variants.vcf.gz; rename to our convention.
  local final_vcf="$out_dir/${prefix}.vcf.gz"
  local final_tbi="$final_vcf.tbi"
  local final_stats="${final_vcf}.stats"

  if ! is_dry_run; then
    if [[ -f "$rundir/results/variants/variants.vcf.gz" ]]; then
      cp "$rundir/results/variants/variants.vcf.gz"     "$final_vcf"
      cp "$rundir/results/variants/variants.vcf.gz.tbi" "$final_tbi"
    fi
  fi

  emit_status "[$sample] Strelka2 germline: bcftools stats"
  module_bcftools_stats --vcf "$final_vcf" --ref "$genome_resolved" --out "$final_stats"

  if ! is_dry_run; then
    [[ -f "$final_vcf" ]]   && emit_artifact variant_vcf       "$final_vcf"
    [[ -f "$final_tbi" ]]   && emit_artifact variant_vcf_index "$final_tbi"
    [[ -f "$final_stats" ]] && emit_artifact variant_stats     "$final_stats"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  outputs_json="$(ngs_state_files_json --kind=variant_vcf "$final_vcf")"
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","exome":"%s","targeted":"%s"}' \
    "$aligner" "${Deduplication:-TRUE}" "${exome:-FALSE}" "${targeted:-FALSE}")"

  ngs_state_stage_end "$sample" Strelka2_germline completed \
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
  pipeline_run Strelka2_germline
fi
