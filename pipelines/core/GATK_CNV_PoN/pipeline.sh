#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/GATK_CNV_PoN/pipeline.sh
#
# Project-scope. Builds a GATK4 CNV Panel of Normals from the cohort's
# normal samples (Group ∈ $pon_groups). Output is consumed by GATK_CNV.
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

_in_csv() {
  # _in_csv "needle" "comma,separated,list"
  local needle="$1" csv="$2"
  local IFS=','
  local item
  for item in $csv; do
    [[ "$item" == "$needle" ]] && return 0
  done
  return 1
}

run_GATK_CNV_PoN_for_project() {
  local aligner="${Aligner:-STAR}"
  local out_dir="$work_dir/$NGS_PROJECT_SENTINEL/GATK_CNV_PoN"
  mkdir -p "$out_dir"
  module_reset

  if [[ ! -f "${cnv_intervals:-}" ]] && ! is_dry_run; then
    log_error "[_project] cnv_intervals not found: ${cnv_intervals:-<unset>}"
    return 1
  fi

  # ---- Pick normal samples from cohort -------------------------------
  local pon_groups="${pon_groups:-Normal}"
  local -a normals=()
  local s
  for s in "${samples[@]}"; do
    local g="${Group_dict[$s]:-}"
    if [[ -z "$g" ]]; then
      log_warn "[_project] $s has no Group; skipping in PoN"
      continue
    fi
    if _in_csv "$g" "$pon_groups"; then
      normals+=("$s")
    fi
  done

  log_info "[_project] PoN candidates (Group ∈ $pon_groups): ${#normals[@]} sample(s)"

  if (( ${#normals[@]} < ${pon_min_samples:-5} )) && ! is_dry_run; then
    log_error "[_project] only ${#normals[@]} normal sample(s) found; need ≥${pon_min_samples:-5}"
    log_error "          (matched Group=$pon_groups in SampleInfoFile)"
    return 1
  fi

  # Record the panel composition for reproducibility.
  local normals_list="$out_dir/normals.txt"
  printf '%s\n' "${normals[@]}" > "$normals_list"

  # ---- Per-normal CollectReadCounts ---------------------------------
  local -a count_files=()
  for s in "${normals[@]}"; do
    local align_dir="$work_dir/$s/Alignment-${aligner}"
    local bam
    bam="$(_resolve_input_bam "$s" "$aligner" "$align_dir")"
    if [[ ! -f "$bam" ]] && ! is_dry_run; then
      log_error "[_project] cannot find normal BAM: $bam"
      return 1
    fi
    local counts="$out_dir/${s}.${aligner}.counts.hdf5"
    emit_status "[_project] PoN: CollectReadCounts $s"
    module_gatk_collect_read_counts \
      --bam "$bam" \
      --intervals "${cnv_intervals:-<intervals>}" \
      --out "$counts"
    count_files+=("$counts")
  done

  # ---- Build PoN ----------------------------------------------------
  local pon="$out_dir/cnv.pon.hdf5"
  emit_status "[_project] PoN: CreateReadCountPanelOfNormals (${#count_files[@]} normals)"
  if (( ${#count_files[@]} > 0 )); then
    module_gatk_create_read_count_pon --out "$pon" \
      $(printf -- '--counts %s ' "${count_files[@]}")
  elif is_dry_run; then
    # Dry-run with no normals: emit a synthetic call for the trace.
    module_gatk_create_read_count_pon --out "$pon" --counts "<normal.counts.hdf5>"
  fi

  if ! is_dry_run; then
    [[ -f "$pon"          ]] && emit_artifact cnv_pon              "$pon"
    [[ -f "$normals_list" ]] && emit_artifact cnv_pon_normals_list "$normals_list"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  if (( ${#count_files[@]} > 0 )); then
    inputs_json="$(ngs_state_files_json --kind=cnv_read_counts "${count_files[@]}")"
  else
    inputs_json='[]'
  fi
  if [[ -f "$pon" ]]; then
    outputs_json="$(ngs_state_files_json --kind=cnv_pon "$pon")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","pon_groups":"%s","n_normals":%d,"intervals":"%s"}' \
    "$aligner" "${Deduplication:-TRUE}" "$pon_groups" "${#normals[@]}" \
    "$(basename "${cnv_intervals:-}")")"

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" GATK_CNV_PoN completed \
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
  pipeline_run GATK_CNV_PoN
fi
