#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/MethylationExtraction/pipeline.sh
#
# Per-sample bismark methylation extraction. The legacy
# GeneralSteps/Alignment.sh ran this inline at the end of bismark
# alignment; the new framework splits the concerns so users can re-run
# extraction with different parameters without re-aligning.
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
module_load bismark

_resolve_input_bam() {
  local sample="$1" aligner="$2" align_dir="$3"
  # Bismark BAM names vary by paired/single + dedup; pick the first match.
  local f
  for pattern in \
    "$align_dir/${sample}.${aligner}.deduplicated.bam" \
    "$align_dir/${sample}.${aligner}.bam" \
    "$align_dir/"*"_pe.deduplicated.bam" \
    "$align_dir/"*"_pe.bam" \
    "$align_dir/"*"_se.deduplicated.bam" \
    "$align_dir/"*"_se.bam"; do
    for f in $pattern; do
      [[ -f "$f" ]] && { echo "$f"; return; }
    done
  done
  # Fallback (dry-run): synthesise the canonical name.
  echo "$align_dir/${sample}.${aligner}.bam"
}

run_MethylationExtraction_for_sample() {
  local sample="${1:?sample required}"
  local aligner="${Aligner:-bismark}"
  local sample_dir="$work_dir/$sample"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local out_dir="$align_dir/MethylationExtraction"
  mkdir -p "$out_dir"
  module_reset

  local bam
  bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ ! -f "$bam" ]] && ! is_dry_run; then
    log_error "[$sample] cannot find Bismark BAM under $align_dir (run Alignment with Aligner=bismark first)"
    return 1
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi
  local genome_folder_resolved="${genome_folder:-${IGenome_Bismark_index:-${IGenome_Genome_dir:-}}}"
  if [[ -z "$genome_folder_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] genome_folder not set and not resolvable from iGenomes"
    return 1
  fi
  [[ -z "$genome_folder_resolved" ]] && genome_folder_resolved="<bismark_genome>"

  local layout="${Layout_dict[$sample]:-PE}"
  local layout_arg
  case "$layout" in
    PE) layout_arg="--paired" ;;
    SE) layout_arg="--single" ;;
    *)  layout_arg="--paired" ;;
  esac

  emit_status "[$sample] MethylationExtraction: bismark_methylation_extractor"
  module_bismark_methylation_extractor \
    --bam "$bam" \
    --out-dir "$out_dir" \
    --genome-folder "$genome_folder_resolved" \
    --threads "${threads:-8}" \
    "$layout_arg"

  # bismark2report — collect the auxiliary reports the legacy pipeline expected.
  local report_dir="$align_dir/bismark2report"
  mkdir -p "$report_dir"

  local align_rpt split_rpt mbias_rpt dedup_rpt
  align_rpt="$(ls "$align_dir"/*_[SP]E_report.txt 2>/dev/null | head -1)"
  split_rpt="$(ls "$out_dir"/*_splitting_report.txt 2>/dev/null | head -1)"
  mbias_rpt="$(ls "$out_dir"/*M-bias.txt          2>/dev/null | head -1)"
  dedup_rpt="$(ls "$align_dir"/*deduplication_report.txt 2>/dev/null | head -1)"

  if [[ -n "$align_rpt" ]] || is_dry_run; then
    emit_status "[$sample] MethylationExtraction: bismark2report"
    module_bismark2report \
      --dir "$report_dir" \
      --alignment-report "${align_rpt:-<alignment_report>}" \
      ${split_rpt:+--splitting-report "$split_rpt"} \
      ${mbias_rpt:+--mbias-report "$mbias_rpt"} \
      ${dedup_rpt:+--dedup-report "$dedup_rpt"}
  else
    log_warn "[$sample] no Bismark _[SP]E_report.txt found in $align_dir — skipping bismark2report"
  fi

  if ! is_dry_run; then
    local cov bg cyto html
    cov="$(ls "$out_dir"/*.bismark.cov.gz       2>/dev/null | head -1)"
    bg="$(ls "$out_dir"/*.bedGraph.gz           2>/dev/null | head -1)"
    cyto="$(ls "$out_dir"/*CpG_report.txt.gz    2>/dev/null | head -1)"
    html="$(ls "$report_dir"/*.html             2>/dev/null | head -1)"
    [[ -n "$cov"  ]] && emit_artifact bismark_coverage         "$cov"
    [[ -n "$bg"   ]] && emit_artifact bismark_bedgraph         "$bg"
    [[ -n "$cyto" ]] && emit_artifact bismark_cytosine_report  "$cyto"
    [[ -n "$html" ]] && emit_artifact bismark_html_report      "$html"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  local cov_path="$out_dir/${sample}.${aligner}.bismark.cov.gz"
  if [[ -f "$cov_path" ]]; then
    outputs_json="$(ngs_state_files_json --kind=bismark_coverage "$cov_path")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"Aligner":"%s","SequenceType":"%s","layout":"%s"}' \
    "$aligner" "${SequenceType:-bisulfite}" "$layout")"

  ngs_state_stage_end "$sample" MethylationExtraction completed \
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
  pipeline_run MethylationExtraction
fi
