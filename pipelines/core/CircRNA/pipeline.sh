#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/CircRNA/pipeline.sh
#
# Per-sample circRNA detection. Runs STAR with chimeric-junction params,
# then CIRCexplorer2 parse + annotate to produce the canonical
# circularRNA_known.txt + back_spliced_junction.bed.
#
# We re-align here (rather than reusing Alignment's BAM) because the chim
# parameters needed for CIRCexplorer2 are different from a standard
# RNA-seq alignment — bundling them keeps CircRNA self-contained.
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
module_load star
module_load circexplorer2

run_CircRNA_for_sample() {
  local sample="${1:?sample required}"
  local sample_dir="$work_dir/$sample"
  local out_dir="$sample_dir/CircRNA"
  local align_dir="$out_dir/star_chim"
  local parse_dir="$out_dir/parse"
  local annot_dir="$out_dir/annotate"
  mkdir -p "$align_dir" "$parse_dir" "$annot_dir"
  module_reset

  local layout="${Layout_dict[$sample]:-PE}"
  local trim_dir="$sample_dir/PreAlignmentQC"
  local fq1 fq2=""
  if [[ "$layout" == "PE" ]]; then
    fq1="$trim_dir/${sample}_trim_1.fq.gz"
    fq2="$trim_dir/${sample}_trim_2.fq.gz"
    if [[ ! -f "$fq1" ]] && ! is_dry_run; then
      log_error "[$sample] expected trimmed FASTQ from preAlignmentQC: $fq1"
      return 1
    fi
  else
    fq1="$trim_dir/${sample}_trim.fq.gz"
  fi

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build" || true
  fi
  local index
  index="$(ngs_aligner_index STAR 2>/dev/null || true)"
  if [[ -z "$index" ]]; then
    if is_dry_run; then
      index="<index_for_STAR>"
    else
      log_error "[$sample] STAR index not found; cannot detect circRNAs"
      return 1
    fi
  fi
  local genome_resolved="${genome:-${IGenome_Genome:-${IGenome_Fasta:-}}}"
  if [[ -z "$genome_resolved" ]] && ! is_dry_run; then
    log_error "[$sample] genome FASTA not set and not resolvable from iGenomes"
    return 1
  fi
  [[ -z "$genome_resolved" ]] && genome_resolved="<genome>"

  if [[ ! -f "${circexplorer2_ref:-}" ]] && ! is_dry_run; then
    log_error "[$sample] circexplorer2_ref not found: ${circexplorer2_ref:-<unset>}"
    return 1
  fi

  # ---- 1. STAR with chim params ------------------------------------
  emit_status "[$sample] CircRNA: STAR (chimSegmentMin=${chim_segment_min:-10})"
  local chim_extra="--chimSegmentMin ${chim_segment_min:-10} --chimOutType Junctions"
  module_star_align \
    --index "$index" \
    --out-prefix "$align_dir/${sample}." \
    --in1 "$fq1" ${fq2:+--in2 "$fq2"} \
    --threads "${threads:-8}" \
    --extra "$chim_extra"

  local chim_junction="$align_dir/${sample}.Chimeric.out.junction"

  # ---- 2. CIRCexplorer2 parse --------------------------------------
  local bsj_bed="$parse_dir/back_spliced_junction.bed"
  emit_status "[$sample] CircRNA: CIRCexplorer2 parse"
  local -a parse_args=(--in "$chim_junction" --out "$bsj_bed" --aligner STAR)
  [[ "$layout" == "PE" ]] && parse_args+=(--paired)
  module_circexplorer2_parse "${parse_args[@]}"

  # ---- 3. CIRCexplorer2 annotate -----------------------------------
  local annot_out="$annot_dir/circularRNA_known.txt"
  emit_status "[$sample] CircRNA: CIRCexplorer2 annotate"
  module_circexplorer2_annotate \
    --bsj "$bsj_bed" \
    --ref "${circexplorer2_ref:-<ref.txt>}" \
    --genome "$genome_resolved" \
    --out "$annot_out"

  if ! is_dry_run; then
    [[ -f "$bsj_bed"  ]] && emit_artifact circrna_bsj_bed "$bsj_bed"
    [[ -f "$annot_out" ]] && emit_artifact circrna_table  "$annot_out"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  if [[ "$layout" == "PE" ]]; then
    inputs_json="$(ngs_state_files_json --kind=trimmed_fastq "$fq1" "$fq2")"
  else
    inputs_json="$(ngs_state_files_json --kind=trimmed_fastq "$fq1")"
  fi
  if [[ -f "$annot_out" ]]; then
    outputs_json="$(ngs_state_files_json --kind=circrna_table "$annot_out")"
  else
    outputs_json='[]'
  fi
  params_json="$(printf '{"chim_segment_min":%s,"layout":"%s","ref":"%s"}' \
    "${chim_segment_min:-10}" "$layout" "$(basename "${circexplorer2_ref:-}")")"

  ngs_state_stage_end "$sample" CircRNA completed \
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
  pipeline_run CircRNA
fi
