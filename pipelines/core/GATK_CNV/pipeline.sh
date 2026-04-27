#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/GATK_CNV/pipeline.sh
#
# GATK4 CNV calling — best-practices germline / tumor-only flow against a
# pre-built Panel of Normals. Replaces the broken legacy GATK CNV scripts
# (which were copy-paste duplicates of the short-variant pipelines).
#
# Per-sample stages:
#   1. CollectReadCounts        bam → hdf5
#   2. DenoiseReadCounts        hdf5 + PoN → standardised + denoised TSVs
#   3. ModelSegments            denoised → segments
#   4. CallCopyRatioSegments    segments → called.seg (+/- amplification, deletion)
#   5. PlotModeledSegments      (optional) PNG visualisation
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

run_GATK_CNV_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-STAR}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/GATK_CNV"
  mkdir -p "$out_dir"
  module_reset

  local bam
  bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ ! -f "$bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find input BAM: $bam (run Alignment first)"
    return 1
  fi

  if [[ ! -f "${cnv_intervals:-}" ]] && ! is_dry_run; then
    log_error "[$sample] cnv_intervals not found: ${cnv_intervals:-<unset>}"
    return 1
  fi
  if [[ ! -f "${cnv_pon:-}" ]] && ! is_dry_run; then
    log_error "[$sample] cnv_pon not found: ${cnv_pon:-<unset>} (build with GATK_CNV_PoN first)"
    return 1
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi
  local genome_resolved="${genome:-${IGenome_Genome:-${IGenome_Fasta:-}}}"
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"

  local prefix
  prefix="$(basename "${bam%.bam}").gatk"
  local counts="$out_dir/${prefix}.counts.hdf5"
  local std_cr="$out_dir/${prefix}.standardizedCR.tsv"
  local denoised_cr="$out_dir/${prefix}.denoisedCR.tsv"
  local segments_dir="$out_dir/segments"
  local called_seg="$out_dir/${prefix}.called.seg"
  mkdir -p "$segments_dir"

  # ---- 1. CollectReadCounts ----
  emit_status "[$sample] GATK CNV: CollectReadCounts"
  module_gatk_collect_read_counts \
    --bam "$bam" \
    --intervals "${cnv_intervals:-<intervals>}" \
    --out "$counts"

  # ---- 2. DenoiseReadCounts ----
  emit_status "[$sample] GATK CNV: DenoiseReadCounts"
  module_gatk_denoise_read_counts \
    --counts "$counts" \
    --pon "${cnv_pon:-<pon.hdf5>}" \
    --standardized-out "$std_cr" \
    --denoised-out "$denoised_cr"

  # ---- 3. ModelSegments ----
  emit_status "[$sample] GATK CNV: ModelSegments"
  module_gatk_model_segments \
    --denoised "$denoised_cr" \
    --out-dir "$segments_dir" \
    --out-prefix "${prefix}"

  # ---- 4. CallCopyRatioSegments ----
  emit_status "[$sample] GATK CNV: CallCopyRatioSegments"
  local cr_seg="$segments_dir/${prefix}.cr.seg"
  module_gatk_call_copy_ratio_segments \
    --in "$cr_seg" \
    --out "$called_seg"

  # ---- 5. PlotModeledSegments (optional) ----
  if [[ "${emit_plots:-true}" == "true" ]]; then
    local seq_dict="${genome_resolved%.fa*}.dict"
    if [[ -f "$seq_dict" ]] || is_dry_run; then
      emit_status "[$sample] GATK CNV: PlotModeledSegments"
      module_gatk_plot_modeled_segments \
        --denoised "$denoised_cr" \
        --segments "$cr_seg" \
        --sequence-dict "${seq_dict:-<dict>}" \
        --out-dir "$out_dir/plots" \
        --out-prefix "$prefix"
    else
      log_warn "[$sample] sequence dictionary not found: $seq_dict — skipping plots"
    fi
  fi

  if ! is_dry_run; then
    [[ -f "$counts"      ]] && emit_artifact cnv_read_counts          "$counts"
    [[ -f "$denoised_cr" ]] && emit_artifact cnv_denoised_copy_ratios "$denoised_cr"
    [[ -f "$called_seg"  ]] && emit_artifact cnv_called_segments     "$called_seg"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  if [[ -f "$called_seg" ]]; then
    outputs_json="$(ngs_state_files_json --kind=cnv_called_segments "$called_seg")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"Aligner":"%s","Deduplication":"%s","emit_plots":%s,"intervals":"%s","pon":"%s"}' \
    "$aligner" "${Deduplication:-TRUE}" "${emit_plots:-true}" \
    "$(basename "${cnv_intervals:-}")" "$(basename "${cnv_pon:-}")")"

  ngs_state_stage_end "$sample" GATK_CNV completed \
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
  pipeline_run GATK_CNV
fi
