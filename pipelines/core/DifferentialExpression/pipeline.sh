#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/DifferentialExpression/pipeline.sh
#
# Project-scope differential-expression pipeline. Wraps an edgeR-based R
# script that the user can swap out via $de_script. Schema requires
# group_a and group_b — the comparison is always logFC = group_b / group_a
# so the labels match what the user expects to see in the volcano plot.
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

run_DifferentialExpression_for_project() {
  module_reset

  # Prefer batch-corrected counts when available, fall back to raw merged counts.
  # Users can force one or the other by setting $use_corrected={true,false}.
  local raw_matrix="$work_dir/$NGS_PROJECT_SENTINEL/Quantification/counts.matrix.tab"
  local corrected_matrix="$work_dir/$NGS_PROJECT_SENTINEL/BatchCorrection/counts.matrix.corrected.tab"
  local matrix matrix_kind
  case "${use_corrected:-auto}" in
    true)
      matrix="$corrected_matrix"; matrix_kind="counts_matrix_corrected"
      ;;
    false)
      matrix="$raw_matrix"; matrix_kind="counts_matrix"
      ;;
    *)
      if [[ -f "$corrected_matrix" ]]; then
        matrix="$corrected_matrix"; matrix_kind="counts_matrix_corrected"
        log_info "[_project] using batch-corrected counts: $corrected_matrix"
      else
        matrix="$raw_matrix"; matrix_kind="counts_matrix"
      fi
      ;;
  esac

  local out_dir="$work_dir/$NGS_PROJECT_SENTINEL/DifferentialExpression"
  mkdir -p "$out_dir"

  if [[ ! -f "$matrix" ]] && ! is_dry_run; then
    log_error "[_project] missing counts matrix: $matrix (run MergeCounts first)"
    return 1
  fi
  if [[ -z "${SampleInfoFile:-}" ]] || { [[ ! -f "${SampleInfoFile:-}" ]] && ! is_dry_run; }; then
    log_error "[_project] SampleInfoFile is required for DifferentialExpression (Group column)"
    return 1
  fi

  local de_script="${de_script:-$_PIPELINE_DIR/scripts/de.R}"
  if [[ ! -f "$de_script" ]] && ! is_dry_run; then
    log_error "[_project] de_script not found: $de_script"
    return 1
  fi

  emit_status "[_project] DifferentialExpression: $group_b vs $group_a (FDR<${max_padj}, |log2FC|>=${min_log2fc})"

  if ! module_rscript_run \
       --script "$de_script" \
       --required "edgeR" \
       --args "$matrix $SampleInfoFile $out_dir ${max_padj:-0.05} ${min_log2fc:-1} ${min_count:-10} $group_a $group_b"; then
    log_error "[_project] DifferentialExpression R script failed"
    return 1
  fi

  if ! is_dry_run; then
    local results_full="$out_dir/de_results_${group_b}_vs_${group_a}.tab"
    local results_sig="$out_dir/de_significant_${group_b}_vs_${group_a}.tab"
    local ma_pdf="$out_dir/ma_${group_b}_vs_${group_a}.pdf"
    local volcano_pdf="$out_dir/volcano_${group_b}_vs_${group_a}.pdf"
    local summary_json="$out_dir/summary.json"
    [[ -f "$results_full"  ]] && emit_artifact de_results_full         "$results_full"
    [[ -f "$results_sig"   ]] && emit_artifact de_results_significant  "$results_sig"
    [[ -f "$summary_json"  ]] && emit_artifact de_summary_json         "$summary_json"
    [[ -f "$ma_pdf"        ]] && emit_artifact de_ma_plot              "$ma_pdf"
    [[ -f "$volcano_pdf"   ]] && emit_artifact de_volcano_plot         "$volcano_pdf"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind="$matrix_kind" "$matrix")"
  local -a out_paths=()
  for f in "$out_dir/de_results_${group_b}_vs_${group_a}.tab" \
           "$out_dir/de_significant_${group_b}_vs_${group_a}.tab" \
           "$out_dir/summary.json"; do
    [[ -f "$f" ]] && out_paths+=("$f")
  done
  if (( ${#out_paths[@]} > 0 )); then
    outputs_json="$(ngs_state_files_json --kind=de_results "${out_paths[@]}")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"group_a":"%s","group_b":"%s","max_padj":%s,"min_log2fc":%s,"min_count":%s,"matrix_kind":"%s"}' \
    "${group_a}" "${group_b}" "${max_padj:-0.05}" "${min_log2fc:-1}" "${min_count:-10}" "$matrix_kind")"

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" DifferentialExpression completed \
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
  pipeline_run DifferentialExpression
fi
