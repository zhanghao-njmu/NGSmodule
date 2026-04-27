#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/postAlignmentQC/pipeline.sh
#
# Quality assessment of aligned BAMs from the Alignment pipeline.
# Composes RSeQC, preseq, mosdepth, goleft modules. Skips RNA-specific
# analyses for DNA / bisulfite data.
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
module_load rseqc
module_load preseq
module_load mosdepth
module_load goleft

# ---------------------------------------------------------------------------
# Resolve a BED file from a GTF on demand, caching alongside the GTF so
# every sample reuses it. Mirrors the legacy genes.bed convention.
# ---------------------------------------------------------------------------
_resolve_genes_bed() {
  local gtf="$1"
  [[ -z "$gtf" ]] && { echo ""; return 0; }
  local bed="${gtf%/*}/genes.bed"
  if [[ ! -f "$bed" ]] && ! is_dry_run; then
    if command -v gtfToGenePred >/dev/null 2>&1 && command -v genePredToBed >/dev/null 2>&1; then
      local tmp="${TMPDIR:-/tmp}/ngs-gtf-$$.genePred"
      run_step "gtfToGenePred + genePredToBed" bash -c \
        "gtfToGenePred '$gtf' '$tmp' && genePredToBed '$tmp' '$bed' && rm -f '$tmp'"
    else
      log_warn "gtfToGenePred/genePredToBed not on PATH — RSeQC modules will fail"
    fi
  fi
  echo "$bed"
}

# ---------------------------------------------------------------------------
# Resolve the input BAM by aligner+deduplication state. Mirrors the legacy
# postAlignmentQC.sh logic.
# ---------------------------------------------------------------------------
_resolve_input_bam() {
  local sample="$1" aligner="$2" align_dir="$3"
  if [[ "${SequenceType:-}" == "bisulfite" ]] || [[ "${SequenceType:-}" == "BSdna" ]]; then
    if [[ "$aligner" =~ ^bismark ]]; then
      ls "$align_dir"/*.bam 2>/dev/null | head -1
      return
    fi
  fi
  case "${Deduplication:-NONE}" in
    FALSE) echo "$align_dir/${sample}.${aligner}.markdup.bam" ;;
    TRUE)  echo "$align_dir/${sample}.${aligner}.dedup.bam" ;;
    *)     echo "$align_dir/${sample}.${aligner}.sorted.bam" ;;
  esac
}

# ---------------------------------------------------------------------------
# Per-sample worker.
# ---------------------------------------------------------------------------
run_postAlignmentQC_for_sample() {
  local sample="${1:?sample required}"
  local sample_dir="$work_dir/$sample"
  local aligner="${Aligner:-STAR}"
  local seqtype="${SequenceType:-rna}"
  local align_dir="$sample_dir/Alignment-${aligner}"
  local qc_dir="$align_dir/postAlignmentQC"
  module_reset

  local bam
  bam="$(_resolve_input_bam "$sample" "$aligner" "$align_dir")"
  if [[ -z "$bam" ]] || { [[ ! -f "$bam" ]] && ! is_dry_run; }; then
    log_error "[$sample] cannot find input BAM: $bam (run Alignment first)"
    return 1
  fi

  local genes_bed
  genes_bed="$(_resolve_genes_bed "${gtf:-}")"

  # ---- RSeQC ------------------------------------------------------------
  local rseqc_dir="$qc_dir/RSeQC"
  mkdir -p "$rseqc_dir"
  emit_status "[$sample] postAlignmentQC: RSeQC"

  module_rseqc_bam_stat \
    --bam "$bam" \
    --out "$rseqc_dir/${sample}.${aligner}.bam_stat.txt"

  if [[ -n "$genes_bed" ]]; then
    module_rseqc_infer_experiment \
      --bam "$bam" --bed "$genes_bed" \
      --out "$rseqc_dir/${sample}.${aligner}.infer_experiment.txt"
    module_rseqc_inner_distance \
      --bam "$bam" --bed "$genes_bed" \
      --out-prefix "$rseqc_dir/${sample}.${aligner}"
    module_rseqc_read_distribution \
      --bam "$bam" --bed "$genes_bed" \
      --out "$rseqc_dir/${sample}.${aligner}.read_distribution.txt"
  fi

  module_rseqc_read_duplication \
    --bam "$bam" \
    --out-prefix "$rseqc_dir/${sample}.${aligner}"
  module_rseqc_read_gc \
    --bam "$bam" \
    --out-prefix "$rseqc_dir/${sample}.${aligner}"

  if [[ "$seqtype" == "rna" ]] && [[ -n "$genes_bed" ]]; then
    emit_status "[$sample] postAlignmentQC: RNA-only RSeQC"
    module_rseqc_genebody_coverage \
      --bam "$bam" --bed "$genes_bed" \
      --out-prefix "$rseqc_dir/${sample}.${aligner}"
    module_rseqc_junction_annotation \
      --bam "$bam" --bed "$genes_bed" \
      --out-prefix "$rseqc_dir/${sample}.${aligner}"
    module_rseqc_junction_saturation \
      --bam "$bam" --bed "$genes_bed" \
      --out-prefix "$rseqc_dir/${sample}.${aligner}"
  fi

  # Bismark BAMs lack the ancillary indexes preseq/mosdepth need; skip.
  local skip_extras=0
  if [[ "$seqtype" == "bisulfite" ]] && [[ "$aligner" =~ ^bismark ]]; then
    skip_extras=1
  fi

  if (( ! skip_extras )); then
    # ---- preseq -------------------------------------------------------
    local preseq_out="$qc_dir/Preseq/${sample}.${aligner}.preseq.txt"
    mkdir -p "$(dirname "$preseq_out")"
    emit_status "[$sample] postAlignmentQC: preseq"
    module_preseq_lc_extrap --bam "$bam" --out "$preseq_out"

    # ---- goleft indexcov ---------------------------------------------
    emit_status "[$sample] postAlignmentQC: goleft indexcov"
    module_goleft_indexcov --bam "$bam" --dir "$qc_dir/goleft"

    # ---- mosdepth ----------------------------------------------------
    local mosdepth_dir="$qc_dir/mosdepth"
    mkdir -p "$mosdepth_dir"
    emit_status "[$sample] postAlignmentQC: mosdepth"
    (
      cd "$mosdepth_dir"
      module_mosdepth_run \
        --bam "$bam" \
        --prefix "${sample}.${aligner}" \
        --threads "${threads:-4}"
    )
  fi

  if ! is_dry_run; then
    emit_artifact rseqc_bam_stat       "$rseqc_dir/${sample}.${aligner}.bam_stat.txt"
    [[ -f "$rseqc_dir/${sample}.${aligner}.infer_experiment.txt" ]] && \
      emit_artifact rseqc_infer_experiment "$rseqc_dir/${sample}.${aligner}.infer_experiment.txt"
    if (( ! skip_extras )); then
      emit_artifact preseq_lc_extrap "$qc_dir/Preseq/${sample}.${aligner}.preseq.txt"
      [[ -f "$qc_dir/mosdepth/${sample}.${aligner}.mosdepth.summary.txt" ]] && \
        emit_artifact mosdepth_summary "$qc_dir/mosdepth/${sample}.${aligner}.mosdepth.summary.txt"
      [[ -f "$qc_dir/goleft/index.html" ]] && \
        emit_artifact goleft_indexcov_html "$qc_dir/goleft/index.html"
    fi
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$bam")"
  local out_files=("$rseqc_dir/${sample}.${aligner}.bam_stat.txt")
  (( ! skip_extras )) && out_files+=("$qc_dir/Preseq/${sample}.${aligner}.preseq.txt")
  outputs_json="$(ngs_state_files_json --kind=postAlignmentQC "${out_files[@]}")"
  params_json="$(printf '{"Aligner":"%s","SequenceType":"%s","Deduplication":"%s"}' \
    "$aligner" "$seqtype" "${Deduplication:-NONE}")"

  ngs_state_stage_end "$sample" postAlignmentQC completed \
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
  pipeline_run postAlignmentQC
fi
