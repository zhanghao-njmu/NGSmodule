#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/postQuantificationQC/pipeline.sh
#
# Project-scope pipeline. Runs a single R script against the merged counts
# matrix from MergeCounts. The bundled scripts/qc.R uses base R only, so it
# works on any Rscript install; users can override with --qc_script to plug
# in their own analysis with Bioconductor packages.
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
module_load rscript

run_postQuantificationQC_for_project() {
  module_reset

  local matrix="$work_dir/$NGS_PROJECT_SENTINEL/Quantification/counts.matrix.tab"
  local out_dir="$work_dir/$NGS_PROJECT_SENTINEL/postQuantificationQC"
  mkdir -p "$out_dir"

  if [[ ! -f "$matrix" ]] && ! is_dry_run; then
    log_error "[_project] missing counts matrix: $matrix (run MergeCounts first)"
    return 1
  fi

  local qc_script="${qc_script:-$_PIPELINE_DIR/scripts/qc.R}"
  if [[ ! -f "$qc_script" ]] && ! is_dry_run; then
    log_error "[_project] qc_script not found: $qc_script"
    return 1
  fi

  emit_status "[_project] postQuantificationQC: running $(basename "$qc_script")"

  if ! module_rscript_run \
       --script "$qc_script" \
       --args "$matrix ${SampleInfoFile:-} $out_dir"; then
    log_error "[_project] postQuantificationQC R script failed"
    return 1
  fi

  if ! is_dry_run; then
    [[ -f "$out_dir/library_sizes.tab"     ]] && emit_artifact postQC_library_sizes "$out_dir/library_sizes.tab"
    [[ -f "$out_dir/sample_correlation.tab" ]] && emit_artifact postQC_correlation   "$out_dir/sample_correlation.tab"
    [[ -f "$out_dir/pca.tab"                ]] && emit_artifact postQC_pca           "$out_dir/pca.tab"
    [[ -f "$out_dir/library_sizes.pdf"      ]] && emit_artifact postQC_library_sizes_pdf "$out_dir/library_sizes.pdf"
    [[ -f "$out_dir/sample_correlation.pdf" ]] && emit_artifact postQC_correlation_pdf   "$out_dir/sample_correlation.pdf"
    [[ -f "$out_dir/pca.pdf"                ]] && emit_artifact postQC_pca_pdf           "$out_dir/pca.pdf"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=counts_matrix "$matrix")"
  local -a out_paths=()
  for f in "$out_dir/library_sizes.tab" "$out_dir/sample_correlation.tab" "$out_dir/pca.tab"; do
    [[ -f "$f" ]] && out_paths+=("$f")
  done
  if (( ${#out_paths[@]} > 0 )); then
    outputs_json="$(ngs_state_files_json --kind=postQC_tables "${out_paths[@]}")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"Aligner":"%s","qc_script":"%s","sample_count":%d}' \
    "${Aligner:-STAR}" "$(basename "$qc_script")" "${#samples[@]}")"

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" postQuantificationQC completed \
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
  pipeline_run postQuantificationQC
fi
