#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/DifferentialMethylation/pipeline.sh
#
# Project-scope DM pipeline. Resolves per-sample bismark.cov.gz files
# from each sample's MethylationExtraction output, then hands them to a
# Welch's t-test R script (default) or a user-provided dm_script.
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

run_DifferentialMethylation_for_project() {
  module_reset

  local aligner="${Aligner:-bismark}"
  local out_dir="$work_dir/$NGS_PROJECT_SENTINEL/DifferentialMethylation"
  mkdir -p "$out_dir"

  # Resolve a glob pointing to all per-sample coverage files.
  local cov_glob="$work_dir/*/Alignment-${aligner}/MethylationExtraction/*.bismark.cov.gz"

  if ! is_dry_run; then
    # Quick presence check — at least one file should exist.
    local found
    # shellcheck disable=SC2086  # globbing is intentional
    found="$(ls $cov_glob 2>/dev/null | head -1)"
    if [[ -z "$found" ]]; then
      log_error "[_project] no .bismark.cov.gz files match: $cov_glob (run MethylationExtraction first)"
      return 1
    fi
  fi

  if [[ -z "${SampleInfoFile:-}" ]] || { [[ ! -f "${SampleInfoFile:-}" ]] && ! is_dry_run; }; then
    log_error "[_project] SampleInfoFile required (Group column)"
    return 1
  fi

  local dm_script="${dm_script:-$_PIPELINE_DIR/scripts/dm.R}"
  if [[ ! -f "$dm_script" ]] && ! is_dry_run; then
    log_error "[_project] dm_script not found: $dm_script"
    return 1
  fi

  emit_status "[_project] DifferentialMethylation: $group_b vs $group_a (qval<${max_qval}, |Δ|>=${min_meth_diff}%)"

  if ! module_rscript_run \
       --script "$dm_script" \
       --args "$cov_glob $SampleInfoFile $out_dir $group_a $group_b ${min_cov:-3} ${max_qval:-0.05} ${min_meth_diff:-25}"; then
    log_error "[_project] DifferentialMethylation R script failed"
    return 1
  fi

  if ! is_dry_run; then
    local results_full="$out_dir/dm_results_${group_b}_vs_${group_a}.tab"
    local results_sig="$out_dir/dm_significant_${group_b}_vs_${group_a}.tab"
    local summary_json="$out_dir/summary.json"
    [[ -f "$results_full" ]] && emit_artifact dm_results_full        "$results_full"
    [[ -f "$results_sig"  ]] && emit_artifact dm_results_significant "$results_sig"
    [[ -f "$summary_json" ]] && emit_artifact dm_summary_json        "$summary_json"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  # Inputs: list every cov file we matched (real run only; dry-run = empty)
  if is_dry_run; then
    inputs_json='[]'
  else
    local -a covs=()
    while IFS= read -r f; do covs+=("$f"); done < <(eval "ls $cov_glob 2>/dev/null")
    if (( ${#covs[@]} > 0 )); then
      inputs_json="$(ngs_state_files_json --kind=bismark_coverage "${covs[@]}")"
    else
      inputs_json='[]'
    fi
  fi
  local -a out_paths=()
  for f in "$out_dir/dm_results_${group_b}_vs_${group_a}.tab" \
           "$out_dir/dm_significant_${group_b}_vs_${group_a}.tab" \
           "$out_dir/summary.json"; do
    [[ -f "$f" ]] && out_paths+=("$f")
  done
  if (( ${#out_paths[@]} > 0 )); then
    outputs_json="$(ngs_state_files_json --kind=dm_results "${out_paths[@]}")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"group_a":"%s","group_b":"%s","min_cov":%s,"max_qval":%s,"min_meth_diff":%s}' \
    "${group_a}" "${group_b}" "${min_cov:-3}" "${max_qval:-0.05}" "${min_meth_diff:-25}")"

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" DifferentialMethylation completed \
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
  pipeline_run DifferentialMethylation
fi
