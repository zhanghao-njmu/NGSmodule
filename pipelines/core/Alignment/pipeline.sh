#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/Alignment/pipeline.sh
#
# Reference alignment of trimmed reads. Demonstrates the framework's
# dependency resolution: this pipeline declares `requires: preAlignmentQC`
# in meta.yml, and the orchestrator transparently runs that prerequisite
# for any sample whose state.json shows it incomplete.
#
# Aligner is chosen at runtime via the $Aligner config variable
# (STAR / bwa-mem2 / hisat2 / bismark). For the framework PoC we keep
# the actual command construction minimal and rely on emit_status +
# state.json for tracking.
###############################################################################

set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

# shellcheck source=../../lib/orchestrator.sh
source "$_NGS_LIB_DIR/orchestrator.sh"

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

  local layout="${Layout_dict[$sample]:-PE}"

  # The framework guarantees this prerequisite ran (or we're in --no-deps
  # mode and the user is on their own).
  local trim_dir="$sample_dir/PreAlignmentQC"
  local fq1 fq2
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

  # Resolve aligner index from iGenomes (no-op if Genome_direct is set).
  if [[ -n "${iGenomes_Dir:-}" ]] && [[ -n "${Species:-}" ]] && [[ -n "${Source:-}" ]] && [[ -n "${Build:-}" ]]; then
    ngs_resolve_igenomes "$Species" "$Source" "$Build"
  fi
  local index
  index="$(ngs_aligner_index "$aligner")"
  if [[ -z "$index" ]] && [[ "${NGS_DRY_RUN:-0}" != "1" ]]; then
    log_error "[$sample] aligner index not found for: $aligner"
    return 1
  fi

  local bam_out="$out_dir/${sample}.${aligner}.sorted.bam"
  local flagstat_out="$out_dir/${sample}.${aligner}.flagstat"
  local align_log="$out_dir/${aligner}.log"

  # Build the command per aligner. The PoC supports STAR + BWA-MEM2 +
  # HISAT2 with minimal configs; production use should extend each branch
  # with the full parameter set the lab needs.
  local cmd_args=()
  case "${aligner,,}" in
    star)
      cmd_args=(
        STAR
        --runThreadN "${threads:-4}"
        --genomeDir "$index"
        --readFilesIn "$fq1" "${fq2:-}"
        --readFilesCommand zcat
        --outFileNamePrefix "$out_dir/${sample}."
        --outSAMtype BAM SortedByCoordinate
      )
      ;;
    bwa|bwa-mem2|bwa_mem2)
      cmd_args=(
        bwa-mem2 mem -t "${threads:-4}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
        "$index" "$fq1" "${fq2:-}"
      )
      ;;
    hisat2)
      cmd_args=(
        hisat2 -p "${threads:-4}" -x "$index"
        -1 "$fq1" -2 "${fq2:-/dev/null}"
        -S /dev/stdout
      )
      ;;
    *)
      log_error "[$sample] unsupported aligner: $aligner"
      return 1
      ;;
  esac

  if is_dry_run; then
    printf '::dry-run::%s\n' "${cmd_args[*]}"
    printf '::dry-run::samtools sort -@ %s -o %s -\n' "${threads:-4}" "$bam_out"
    printf '::dry-run::samtools index %s\n' "$bam_out"
  else
    run_step "[$sample] $aligner align" \
      bash -c "${cmd_args[*]} 2>'$align_log' | samtools sort -@ ${threads:-4} -o '$bam_out' -" \
      || { ngs_check_logfile "$sample" "$aligner" "$align_log" postcheck $? || return 1; }
    run_step "[$sample] samtools index" samtools index "$bam_out"
    run_step "[$sample] samtools flagstat" \
      bash -c "samtools flagstat '$bam_out' > '$flagstat_out'"
    emit_artifact aligned_bam "$bam_out"
    emit_artifact alignment_stats "$flagstat_out"
  fi

  # Record the stage in the sample state with rich metadata.
  local tools_json inputs_json outputs_json params_json
  tools_json="$(ngs_state_tools_json "${aligner,,}" samtools)"
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
