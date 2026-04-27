#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/preAlignmentQC/pipeline.sh
#
# Framework-style implementation of the preAlignmentQC step.
#
# This is the PoC for the pipeline framework. Each step under
# pipelines/core/<name>/ provides:
#
#   - meta.yml      → metadata consumed by CLI + backend
#   - pipeline.sh   → defines run_<name>_for_sample() and calls pipeline_run
#   - env.yml       → conda env spec (optional; replaces CheckENV.sh)
#
# Compared to the legacy GeneralSteps/preAlignmentQC.sh (445 lines), the
# orchestration scaffolding is gone — we only describe what a single
# sample's QC looks like. The framework handles concurrency, resume, and
# progress reporting.
#
# Tools used:
#   - fastp (adapter trim, quality filter, paired-read merge)
#   - FastQC (per-base quality + GC content)
###############################################################################

set -euo pipefail

# Locate the framework's lib/ dir relative to this file. Works whether the
# script is invoked directly or sourced by the legacy entry point.
_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

# shellcheck source=../../lib/orchestrator.sh
source "$_NGS_LIB_DIR/orchestrator.sh"

# ---------------------------------------------------------------------------
# Per-sample worker.
#
# Inputs (resolved via ngs_sample_inputs):
#   $work_dir/$sample/run*_${sample}{,_1,_2}.fq.gz
#
# Outputs (under $work_dir/$sample/PreAlignmentQC/):
#   ${sample}_trim_1.fq.gz / ${sample}_trim_2.fq.gz   (PE)
#   ${sample}_trim.fq.gz                              (SE)
#   ${sample}_fastp.html                              (interactive QC report)
#   ${sample}_fastp.json                              (machine-readable metrics)
#   fastqc/                                            (FastQC reports per file)
#
# Idempotent: rerun-safe via ngs_check_logfile.
# ---------------------------------------------------------------------------
run_preAlignmentQC_for_sample() {
  local sample="${1:?sample required}"
  local sample_dir="$work_dir/$sample"
  local out_dir="$sample_dir/PreAlignmentQC"
  mkdir -p "$out_dir"

  local layout="${Layout_dict[$sample]:-PE}"
  local r1_files=() r2_files=()
  if [[ "$layout" == "PE" ]]; then
    ngs_sample_inputs "$sample_dir" PE r1_files r2_files
  else
    ngs_sample_inputs "$sample_dir" SE r1_files
  fi

  if (( ${#r1_files[@]} == 0 )); then
    if [[ "${NGS_DRY_RUN:-0}" == "1" ]]; then
      # Synthetic placeholder so dry-run / `ngsmodule test` works on bare CI.
      r1_files=("$sample_dir/<placeholder>_R1.fq.gz")
      [[ "$layout" == "PE" ]] && r2_files=("$sample_dir/<placeholder>_R2.fq.gz")
    else
      log_error "[$sample] no input FASTQ found under $sample_dir"
      return 1
    fi
  fi

  emit_status "[$sample] preAlignmentQC: validating input"

  # Concatenate multi-run files via process substitution (no extra disk).
  local fq1_in="$out_dir/${sample}.r1.fq.gz"
  local fq2_in="$out_dir/${sample}.r2.fq.gz"

  # We could pipe directly into fastp, but disk-cached intermediates make
  # the resume story simpler — if fastp fails halfway we still have the
  # concatenated input.
  if [[ ! -f "$fq1_in" ]] || [[ "$force" == "TRUE" ]]; then
    log_info "[$sample] concatenating ${#r1_files[@]} R1 file(s)"
    cat "${r1_files[@]}" > "$fq1_in"
    if [[ "$layout" == "PE" ]]; then
      cat "${r2_files[@]}" > "$fq2_in"
    fi
  fi

  emit_status "[$sample] preAlignmentQC: running fastp"
  local fastp_log="$out_dir/fastp.log"
  local fastp_html="$out_dir/${sample}_fastp.html"
  local fastp_json="$out_dir/${sample}_fastp.json"
  local fastp_out1="$out_dir/${sample}_trim_1.fq.gz"
  local fastp_out2="$out_dir/${sample}_trim_2.fq.gz"
  local fastp_out_se="$out_dir/${sample}_trim.fq.gz"

  if ngs_check_logfile "$sample" fastp "$fastp_log" precheck; then
    log_info "[$sample] fastp already done — skip"
  else
    local fastp_args=(
      --thread "${threads_fastp:-${threads:-4}}"
      --json "$fastp_json"
      --html "$fastp_html"
      --qualified_quality_phred "${min_quality:-20}"
      --length_required "${min_length:-50}"
    )
    if [[ "$layout" == "PE" ]]; then
      fastp_args+=(
        --in1 "$fq1_in" --in2 "$fq2_in"
        --out1 "$fastp_out1" --out2 "$fastp_out2"
        --detect_adapter_for_pe
      )
    else
      fastp_args+=( --in1 "$fq1_in" --out1 "$fastp_out_se" )
    fi

    if is_dry_run; then
      printf '::dry-run::fastp %s\n' "${fastp_args[*]}"
      printf 'NGSmodule finished the job [fastp]\n' > "$fastp_log"
    else
      run_step "[$sample] fastp" fastp "${fastp_args[@]}" >"$fastp_log" 2>&1
      ngs_check_logfile "$sample" fastp "$fastp_log" postcheck $? || return 1
    fi
    emit_artifact qc_report_html "$fastp_html"
    emit_artifact qc_report_json "$fastp_json"
  fi

  emit_status "[$sample] preAlignmentQC: running FastQC"
  local fastqc_dir="$out_dir/fastqc"
  local fastqc_log="$out_dir/fastqc.log"
  if ngs_check_logfile "$sample" fastqc "$fastqc_log" precheck; then
    log_info "[$sample] fastqc already done — skip"
  else
    mkdir -p "$fastqc_dir"
    local fastqc_inputs=()
    if [[ "$layout" == "PE" ]]; then
      fastqc_inputs=("$fastp_out1" "$fastp_out2")
    else
      fastqc_inputs=("$fastp_out_se")
    fi
    if is_dry_run; then
      printf '::dry-run::fastqc -t %s -o %s %s\n' \
        "${threads:-4}" "$fastqc_dir" "${fastqc_inputs[*]}"
      printf 'NGSmodule finished the job [fastqc]\n' > "$fastqc_log"
    else
      run_step "[$sample] fastqc" fastqc -t "${threads:-4}" -o "$fastqc_dir" \
        "${fastqc_inputs[@]}" >"$fastqc_log" 2>&1
      ngs_check_logfile "$sample" fastqc "$fastqc_log" postcheck $? || return 1
    fi
    emit_artifact qc_report_html "$fastqc_dir"
  fi

  # Optional: extract a couple of headline metrics for the UI.
  if [[ -f "$fastp_json" ]] && command -v jq >/dev/null 2>&1; then
    local total_reads dup_rate q30
    total_reads=$(jq -r '.summary.before_filtering.total_reads // 0' "$fastp_json")
    dup_rate=$(jq -r '.duplication.rate // 0' "$fastp_json")
    q30=$(jq -r '.summary.before_filtering.q30_rate // 0' "$fastp_json")
    emit_metric "${sample}.total_reads" "$total_reads"
    emit_metric "${sample}.duplication_rate" "$dup_rate"
    emit_metric "${sample}.q30_rate" "$q30"
  fi

  return 0
}

# Allow this file to be run directly: `bash pipelines/core/preAlignmentQC/pipeline.sh`
# When sourced by the legacy entry, the caller dispatches `pipeline_run`
# itself, so we only auto-run when this is the entry point.
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
  pipeline_run preAlignmentQC
fi
