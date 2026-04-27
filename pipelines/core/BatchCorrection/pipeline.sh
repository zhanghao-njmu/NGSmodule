#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/BatchCorrection/pipeline.sh
#
# Project-scope batch correction. Inserts cleanly between MergeCounts and
# DifferentialExpression: read the merged counts → run sva::ComBat_seq (or
# limma::removeBatchEffect) → write a corrected matrix in the same schema.
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

run_BatchCorrection_for_project() {
  module_reset

  local matrix="$work_dir/$NGS_PROJECT_SENTINEL/Quantification/counts.matrix.tab"
  local out_dir="$work_dir/$NGS_PROJECT_SENTINEL/BatchCorrection"
  mkdir -p "$out_dir"

  if [[ ! -f "$matrix" ]] && ! is_dry_run; then
    log_error "[_project] missing counts matrix: $matrix (run MergeCounts first)"
    return 1
  fi
  if [[ -z "${SampleInfoFile:-}" ]] || { [[ ! -f "${SampleInfoFile:-}" ]] && ! is_dry_run; }; then
    log_error "[_project] SampleInfoFile required for BatchCorrection (Batch column)"
    return 1
  fi

  local bc_script="${bc_script:-$_PIPELINE_DIR/scripts/bc.R}"
  if [[ ! -f "$bc_script" ]] && ! is_dry_run; then
    log_error "[_project] bc_script not found: $bc_script"
    return 1
  fi

  local method="${method:-combat_seq}"
  local required_pkg
  case "$method" in
    combat_seq) required_pkg="sva" ;;
    limma)      required_pkg="limma" ;;
    *) log_error "[_project] unknown method=$method (valid: combat_seq, limma)"; return 1 ;;
  esac

  emit_status "[_project] BatchCorrection: method=$method"

  module_rscript_run \
    --script "$bc_script" \
    --required "$required_pkg" \
    --args "$matrix $SampleInfoFile $out_dir $method"

  if ! is_dry_run; then
    [[ -f "$out_dir/counts.matrix.corrected.tab" ]] && \
      emit_artifact counts_matrix_corrected "$out_dir/counts.matrix.corrected.tab"
    [[ -f "$out_dir/pca_before.pdf" ]] && emit_artifact bc_pca_before    "$out_dir/pca_before.pdf"
    [[ -f "$out_dir/pca_after.pdf"  ]] && emit_artifact bc_pca_after     "$out_dir/pca_after.pdf"
    [[ -f "$out_dir/summary.json"   ]] && emit_artifact bc_summary_json  "$out_dir/summary.json"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=counts_matrix "$matrix")"
  if [[ -f "$out_dir/counts.matrix.corrected.tab" ]]; then
    outputs_json="$(ngs_state_files_json --kind=counts_matrix_corrected "$out_dir/counts.matrix.corrected.tab")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"Aligner":"%s","method":"%s","sample_count":%d}' \
    "${Aligner:-STAR}" "$method" "${#samples[@]}")"

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" BatchCorrection completed \
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
  pipeline_run BatchCorrection
fi
