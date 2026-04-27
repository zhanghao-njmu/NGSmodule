#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/checkpoint.sh
#
# Logfile-based resume / completion detection.
#
# Why logfile-based instead of content-hash (Nextflow) or file-presence
# (Snakemake)?
#   - Catches tools that exit 0 but leave corrupt output (we grep for
#     EOFException, "no such file", etc. in their stderr).
#   - Doesn't false-positive when an upstream FASTQ is touched / re-mounted
#     (Nextflow `-resume` blows the entire DAG when input mtime changes).
#   - Trivially extensible: tools add their own error patterns without
#     changing the framework.
#
# Two-tier API:
#   ngs_check_logfile <sample> <tool> <logfile> <mode> [exit_code]
#     mode = precheck   → return 0 (skip) | 1 (run)
#     mode = postcheck  → return 0 (ok)   | 1 (failed)
#
#   ngs_globalcheck_logfile <dir> <logfile_array_var> <force> <sample>
#     Sweeps a sample's working dir and removes any logfile that is
#     "uncompleted with errors" so the next precheck retries it.
#     `force=TRUE` drops every logfile (full rebuild).
#
# Patterns are extensible. Set NGS_ERROR_PATTERN / NGS_COMPLETE_PATTERN
# from a step script to add tool-specific detections.
###############################################################################

# Default patterns. Step scripts may override before sourcing this file
# *or* by reassigning the variable inline.
: "${NGS_ERROR_PATTERN:=(error)|(fatal)|(corrupt)|(interrupt)|(EOFException)|(no such file or directory)}"
: "${NGS_COMPLETE_PATTERN:=(NGSmodule finished the job)}"

# Internal: extract distinct error/complete matches from a logfile.
_ngs_match_log() {
  local logfile="$1" pattern="$2"
  [[ -f "$logfile" ]] || { echo ""; return; }
  grep -ioP "$pattern" "$logfile" 2>/dev/null | sort -u | paste -sd "|"
}

# Append the standard completion marker. Step scripts call this on success.
ngs_mark_complete() {
  local logfile="$1" tool="${2:-step}"
  printf 'NGSmodule finished the job [%s]\n' "$tool" >> "$logfile"
}

# ---------------------------------------------------------------------------
# ngs_check_logfile <sample> <tool> <logfile> <mode> [exit_status]
#
# precheck:
#   - logfile not present              → 1 (start fresh)
#   - logfile has completion marker    → 0 (skip)
#   - logfile has errors               → 1 (rerun)
#   - logfile exists but ambiguous     → 1 (rerun, safe default)
#
# postcheck:
#   - non-zero exit_status             → 1 (failed)
#   - completion marker found          → 0 (idempotent skip)
#   - error pattern in log             → 1 (failed)
#   - no errors and tool exited 0      → 0 (mark complete + return success)
# ---------------------------------------------------------------------------
ngs_check_logfile() {
  local sample="$1" tool="$2" logfile="$3" mode="$4" status="${5:-}"
  local error complete

  if [[ "$mode" == "postcheck" ]] && [[ -n "$status" ]] && [[ "$status" != "0" ]]; then
    return 1
  fi

  if [[ ! -f "$logfile" ]]; then
    [[ "$mode" == "precheck" ]] && return 1
    return 1   # postcheck: missing log = problem
  fi

  complete=$(_ngs_match_log "$logfile" "$NGS_COMPLETE_PATTERN")
  error=$(_ngs_match_log "$logfile" "$NGS_ERROR_PATTERN")

  if [[ -n "$complete" ]]; then
    if [[ "$mode" == "postcheck" ]]; then
      ngs_mark_complete "$logfile" "$tool"
    fi
    return 0
  fi

  if [[ -n "$error" ]]; then
    return 1
  fi

  # No marker, no errors. precheck → safer to rerun. postcheck → tool
  # exited 0 with no detectable error, treat as success and mark.
  if [[ "$mode" == "postcheck" ]]; then
    ngs_mark_complete "$logfile" "$tool"
    return 0
  fi
  return 1
}

# ---------------------------------------------------------------------------
# Sweep a working directory for stale logs from a previous run.
#
# Args:
#   $1  workdir to scan
#   $2  bash array variable name (e.g. logfiles=(fastp.log fastqc.log); pass logfiles[@])
#   $3  force: "TRUE" to delete ALL matching logs, anything else = only
#       delete logs that have errors AND no completion marker
#   $4  sample id (for log messages)
#
# Use this before invoking a step on a sample so the precheck logic sees
# a clean slate where appropriate.
# ---------------------------------------------------------------------------
ngs_globalcheck_logfile() {
  local workdir="$1"
  local -n _names="$2"   # nameref → array of basenames
  local force="$3"
  local sample="$4"

  local find_args=()
  local first=1
  for name in "${_names[@]}"; do
    if (( first )); then
      find_args+=( -name "$name" )
      first=0
    else
      find_args+=( -o -name "$name" )
    fi
  done
  (( ${#find_args[@]} == 0 )) && return 0

  while IFS= read -r log; do
    if [[ "$force" == "TRUE" ]]; then
      rm -f "$log"
      continue
    fi
    if grep -iqP "$NGS_ERROR_PATTERN" "$log" 2>/dev/null && \
       ! grep -iqP "$NGS_COMPLETE_PATTERN" "$log" 2>/dev/null; then
      rm -f "$log"
    fi
  done < <(find "$workdir" \( "${find_args[@]}" \) -type f 2>/dev/null)
}

# Convenience: a pipeline-level "is this sample fully done?" check that
# orchestrators use to decide whether to enqueue the sample at all.
#
# A pipeline is considered complete for a sample when a top-level marker
# file exists. The orchestrator writes this when all its steps succeed.
ngs_pipeline_completed() {
  local sample_workdir="$1" pipeline="$2"
  [[ -f "${sample_workdir}/.ngs_complete.${pipeline}" ]]
}

ngs_pipeline_mark_complete() {
  local sample_workdir="$1" pipeline="$2"
  mkdir -p "$sample_workdir"
  printf 'completed_at=%s\npipeline=%s\nhost=%s\n' \
    "$(date -Iseconds)" "$pipeline" "$(hostname 2>/dev/null || echo unknown)" \
    > "${sample_workdir}/.ngs_complete.${pipeline}"
}

ngs_pipeline_clear_complete() {
  local sample_workdir="$1" pipeline="$2"
  rm -f "${sample_workdir}/.ngs_complete.${pipeline}"
}
