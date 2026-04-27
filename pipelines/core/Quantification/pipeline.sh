#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/Quantification/pipeline.sh
#
# Per-sample read counting via featureCounts. Outputs a counts table and a
# summary file alongside the aligned BAM. Project-level matrix integration
# is handled by a downstream analysis pipeline (not migrated yet).
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
module_load featurecounts

# ---------------------------------------------------------------------------
# Per-sample worker.
#
# Inputs:
#   $work_dir/$sample/Alignment-${Aligner}/${sample}.${Aligner}.sorted.bam
#
# Outputs:
#   $work_dir/$sample/Alignment-${Aligner}/Quantification/
#     ${sample}.${Aligner}.counts.tab
#     ${sample}.${Aligner}.counts.tab.summary
# ---------------------------------------------------------------------------
run_Quantification_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local quant_dir="$align_dir/Quantification"
  module_reset

  local bam="$align_dir/${sample}.${aligner}.sorted.bam"
  if [[ ! -f "$bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find input BAM: $bam (run Alignment first)"
    return 1
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi

  local gtf_resolved="${gtf:-${IGenome_GTF:-}}"
  if [[ -z "$gtf_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] gtf not set and not resolvable from iGenomes"
    return 1
  fi
  if [[ ! -f "$gtf_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] gtf file not found: $gtf_resolved"
    return 1
  fi

  mkdir -p "$quant_dir"
  emit_status "[$sample] Quantification: featureCounts"

  local layout="${Layout_dict[$sample]:-PE}"
  local counts_out="$quant_dir/${sample}.${aligner}.counts.tab"

  local args=(
    --gtf "${gtf_resolved:-<gtf>}"
    --out "$counts_out"
    --bam "$bam"
    --threads "${threads_featurecounts:-${threads:-4}}"
    --strand "${strand_specific:-0}"
    --feature-type "${feature_type:-exon}"
    --meta-feature "${meta_feature:-gene_id}"
  )
  [[ "$layout" == "PE" ]] && args+=(--paired)

  module_featurecounts_count "${args[@]}"

  if ! is_dry_run; then
    emit_artifact counts_tab "$counts_out"
    [[ -f "${counts_out}.summary" ]] && emit_artifact counts_summary "${counts_out}.summary"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  outputs_json="$(ngs_state_files_json --kind=counts_tab "$counts_out")"
  params_json="$(printf '{"Aligner":"%s","strand_specific":"%s","feature_type":"%s","meta_feature":"%s"}' \
    "$aligner" "${strand_specific:-0}" "${feature_type:-exon}" "${meta_feature:-gene_id}")"

  ngs_state_stage_end "$sample" Quantification completed \
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
  pipeline_run Quantification
fi
