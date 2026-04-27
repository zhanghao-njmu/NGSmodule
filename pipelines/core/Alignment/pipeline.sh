#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/Alignment/pipeline.sh
#
# Reference alignment of trimmed reads. Demonstrates two framework
# capabilities at once:
#
#   1. Dependency resolution — declares `requires: preAlignmentQC` in
#      meta.yml; the orchestrator transparently runs that prerequisite
#      for any sample whose state.json shows it incomplete.
#
#   2. Module composition — the per-aligner branch is a few lines of
#      module_<tool>_<action> calls instead of inline command construction.
#      See pipelines/modules/{star,bwa_mem2,hisat2,samtools}.sh.
###############################################################################

set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

# shellcheck source=../../lib/orchestrator.sh
source "$_NGS_LIB_DIR/orchestrator.sh"

# Pipelines compose tool-level modules instead of inlining commands.
# See pipelines/modules/_lib.sh for the contract.
_NGS_MODULES_DIR="$(cd "$_PIPELINE_DIR/../../modules" && pwd)"
# shellcheck source=../../modules/_lib.sh
source "$_NGS_MODULES_DIR/_lib.sh"
module_load star
module_load bwa_mem2
module_load hisat2
module_load samtools

# ---------------------------------------------------------------------------
# Per-sample worker.
#
# Inputs (resolved from prior preAlignmentQC stage):
#   $work_dir/$sample/PreAlignmentQC/${sample}_trim_{1,2}.fq.gz   (PE)
#   $work_dir/$sample/PreAlignmentQC/${sample}_trim.fq.gz         (SE)
#
# Outputs:
#   $work_dir/$sample/Alignment-${Aligner}/${sample}.${Aligner}.sorted.bam
# ---------------------------------------------------------------------------
run_Alignment_for_sample() {
  local sample="${1:?sample required}"
  local sample_dir="$work_dir/$sample"
  local aligner="${Aligner:-STAR}"
  local out_dir="$sample_dir/Alignment-${aligner}"
  mkdir -p "$out_dir"
  module_reset  # fresh per-sample tool tally

  local layout="${Layout_dict[$sample]:-PE}"

  local trim_dir="$sample_dir/PreAlignmentQC"
  local fq1 fq2=""
  if [[ "$layout" == "PE" ]]; then
    fq1="$trim_dir/${sample}_trim_1.fq.gz"
    fq2="$trim_dir/${sample}_trim_2.fq.gz"
    if [[ ! -f "$fq1" ]] && [[ "${NGS_DRY_RUN:-0}" != "1" ]]; then
      log_error "[$sample] expected trimmed FASTQ from preAlignmentQC: $fq1"
      return 1
    fi
  else
    fq1="$trim_dir/${sample}_trim.fq.gz"
  fi

  emit_status "[$sample] Alignment ($aligner): mapping reads"

  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build"
  fi
  local index
  index="$(ngs_aligner_index "$aligner")"
  if [[ -z "$index" ]]; then
    if [[ "${NGS_DRY_RUN:-0}" == "1" ]]; then
      index="<index_for_${aligner}>"   # placeholder so module arg-checks pass
    else
      log_error "[$sample] aligner index not found for: $aligner"
      return 1
    fi
  fi

  local bam_out="$out_dir/${sample}.${aligner}.sorted.bam"
  local flagstat_out="$out_dir/${sample}.${aligner}.flagstat"

  case "${aligner,,}" in
    star)
      module_star_align \
        --index "$index" \
        --out-prefix "$out_dir/${sample}." \
        --in1 "$fq1" ${fq2:+--in2 "$fq2"} \
        --threads "${threads:-4}"
      # STAR writes ${prefix}Aligned.sortedByCoord.out.bam — link to canonical name.
      if ! is_dry_run; then
        local star_bam="$out_dir/${sample}.Aligned.sortedByCoord.out.bam"
        [[ -f "$star_bam" ]] && ln -sf "$(basename "$star_bam")" "$bam_out"
      fi
      ;;
    bwa|bwa-mem2|bwa_mem2)
      local sam="$out_dir/${sample}.bwa.sam"
      if is_dry_run; then
        module_bwa_mem2_align --index "$index" --in1 "$fq1" ${fq2:+--in2 "$fq2"} \
          --rg-sample "$sample" --threads "${threads:-4}"
        module_samtools_sort --in "$sam" --out "$bam_out" --threads "${threads:-4}"
      else
        # Pipe align → sort, no intermediate SAM.
        run_step "[$sample] bwa-mem2 | sort" bash -c "
          bwa-mem2 mem -t '${threads:-4}' \
            -R '@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA' \
            '$index' '$fq1' ${fq2:+'$fq2'} \
          | samtools sort -@ '${threads:-4}' -o '$bam_out' -
        "
        _MODULES_TOOLS_USED[bwa-mem2]=1
        _MODULES_TOOLS_USED[samtools]=1
      fi
      ;;
    hisat2)
      local sam="$out_dir/${sample}.hisat2.sam"
      module_hisat2_align --index "$index" --in1 "$fq1" ${fq2:+--in2 "$fq2"} \
        --out-sam "$sam" --threads "${threads:-4}"
      module_samtools_sort --in "$sam" --out "$bam_out" --threads "${threads:-4}"
      ! is_dry_run && rm -f "$sam"
      ;;
    *)
      log_error "[$sample] unsupported aligner: $aligner"
      return 1
      ;;
  esac

  module_samtools_index --bam "$bam_out"
  module_samtools_flagstat --bam "$bam_out" --out "$flagstat_out"

  if ! is_dry_run; then
    emit_artifact aligned_bam "$bam_out"
    emit_artifact alignment_stats "$flagstat_out"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  if [[ "$layout" == "PE" ]]; then
    inputs_json="$(ngs_state_files_json --kind=trimmed_fastq "$fq1" "$fq2")"
  else
    inputs_json="$(ngs_state_files_json --kind=trimmed_fastq "$fq1")"
  fi
  outputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam_out")"
  params_json="$(printf '{"Aligner":"%s","SequenceType":"%s","threads":%s}' \
    "$aligner" "${SequenceType:-rna}" "${threads:-4}")"

  ngs_state_stage_end "$sample" Alignment completed \
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
  pipeline_run Alignment
fi
