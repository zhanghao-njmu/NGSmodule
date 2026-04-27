#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/core/MergeCounts/pipeline.sh
#
# Project-scope pipeline (scope:project in meta.yml). Runs once across the
# cohort, aggregating per-sample featureCounts output into a single matrix.
#
# featureCounts output schema (from pipelines/modules/featurecounts.sh):
#   # Program:featureCounts ...
#   Geneid<TAB>Chr<TAB>Start<TAB>End<TAB>Strand<TAB>Length<TAB><BAM>
#   <gene1><TAB>...<TAB><count>
#
# Merge strategy: read column 1 (Geneid) + last column (count) per sample,
# join on Geneid using awk's associative arrays. First sample sets the gene
# order; later samples are looked up.
###############################################################################

set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

# shellcheck source=../../lib/orchestrator.sh
source "$_NGS_LIB_DIR/orchestrator.sh"

run_MergeCounts_for_project() {
  local aligner="${Aligner:-STAR}"
  local proj_dir="$work_dir/$NGS_PROJECT_SENTINEL/Quantification"
  mkdir -p "$proj_dir"
  local matrix="$proj_dir/counts.matrix.tab"
  local summary="$proj_dir/counts.summary.merged.tab"

  # ---- Locate per-sample counts files --------------------------------
  local -a count_files=() summary_files=() found_samples=()
  local s f sf
  for s in "${samples[@]}"; do
    f="$work_dir/$s/Alignment-${aligner}/Quantification/${s}.${aligner}.counts.tab"
    sf="${f}.summary"
    if [[ -f "$f" ]] || is_dry_run; then
      count_files+=("$f")
      found_samples+=("$s")
    else
      log_warn "[_project] missing counts.tab for sample: $s ($f)"
    fi
    [[ -f "$sf" ]] && summary_files+=("$sf")
  done

  if (( ${#count_files[@]} == 0 )) && ! is_dry_run; then
    log_error "[_project] no counts files found under $work_dir/<sample>/Alignment-${aligner}/Quantification/"
    return 1
  fi

  emit_status "[_project] MergeCounts: merging ${#count_files[@]} per-sample tables"

  # ---- Merge counts --------------------------------------------------
  if is_dry_run; then
    printf '::dry-run::merge %d counts files → %s\n' "${#count_files[@]}" "$matrix"
  else
    # awk script: build a map Geneid → "count1<TAB>count2<TAB>..." in
    # first-file gene order. Last column of each input row is the count.
    awk -v samples_csv="$(IFS=,; echo "${found_samples[*]}")" '
      BEGIN { FS = OFS = "\t" }
      FNR == 1 {
        file_idx++
        if (file_idx == 1) {
          # First file establishes gene order.
          first_file = 1
        } else {
          first_file = 0
        }
        next
      }
      /^#/ { next }
      first_file && FNR == 2 {
        # featureCounts header: Geneid Chr Start End Strand Length <BAM>
        # We just skip it (already past the comment).
        next
      }
      first_file {
        # Remember gene order + per-gene metadata (Chr,Start,End,Strand,Length).
        gene = $1
        if (!(gene in seen)) {
          gene_order[++n_genes] = gene
          seen[gene] = 1
          # Build metadata from columns 2..6 (Chr,Start,End,Strand,Length).
          meta[gene] = $2 OFS $3 OFS $4 OFS $5 OFS $6
        }
        counts[gene, 1] = $NF
        next
      }
      {
        gene = $1
        counts[gene, file_idx] = $NF
      }
      END {
        n_samples = file_idx
        # Header: Geneid Chr Start End Strand Length sample1 sample2 ...
        split(samples_csv, names, ",")
        printf "Geneid" OFS "Chr" OFS "Start" OFS "End" OFS "Strand" OFS "Length"
        for (i = 1; i <= n_samples; i++) printf OFS "%s", names[i]
        printf "\n"
        for (i = 1; i <= n_genes; i++) {
          gene = gene_order[i]
          printf "%s" OFS "%s", gene, meta[gene]
          for (j = 1; j <= n_samples; j++) {
            v = counts[gene, j]
            printf OFS "%s", (v == "" ? "0" : v)
          }
          printf "\n"
        }
      }
    ' "${count_files[@]}" > "$matrix"

    log_ok "[_project] wrote merged matrix: $matrix"
    emit_artifact counts_matrix "$matrix"
  fi

  # ---- Merge summaries (also two-column data: Status + count) -------
  if (( ${#summary_files[@]} > 0 )) && ! is_dry_run; then
    awk -v samples_csv="$(IFS=,; echo "${found_samples[*]}")" '
      BEGIN { FS = OFS = "\t" }
      FNR == 1 {
        file_idx++
        if (file_idx == 1) first_file = 1; else first_file = 0
        next
      }
      first_file {
        status = $1
        if (!(status in seen)) {
          status_order[++n_statuses] = status
          seen[status] = 1
        }
        vals[status, 1] = $2
        next
      }
      {
        vals[$1, file_idx] = $2
      }
      END {
        n_samples = file_idx
        split(samples_csv, names, ",")
        printf "Status"
        for (i = 1; i <= n_samples; i++) printf OFS "%s", names[i]
        printf "\n"
        for (i = 1; i <= n_statuses; i++) {
          status = status_order[i]
          printf "%s", status
          for (j = 1; j <= n_samples; j++) {
            v = vals[status, j]
            printf OFS "%s", (v == "" ? "0" : v)
          }
          printf "\n"
        }
      }
    ' "${summary_files[@]}" > "$summary"
    log_ok "[_project] wrote merged summary: $summary"
    emit_artifact counts_summary_merged "$summary"
  fi

  # ---- Provenance: record what we merged into state.json ------------
  local inputs_json outputs_json params_json tools_json
  if (( ${#count_files[@]} > 0 )); then
    inputs_json="$(ngs_state_files_json --kind=counts_tab "${count_files[@]}")"
  else
    inputs_json='[]'
  fi
  outputs_json="$(ngs_state_files_json --kind=counts_matrix "$matrix")"
  params_json="$(printf '{"Aligner":"%s","sample_count":%d}' \
    "$aligner" "${#found_samples[@]}")"
  tools_json='[{"name":"awk"}]'

  ngs_state_stage_end "$NGS_PROJECT_SENTINEL" MergeCounts completed \
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
  pipeline_run MergeCounts
fi
