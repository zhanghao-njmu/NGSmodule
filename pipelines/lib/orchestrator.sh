#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/orchestrator.sh
#
# The framework's heart. Pipelines under pipelines/core/ implement a single
# function `run_<pipeline>_for_sample <sample>` and the orchestrator handles
# everything else:
#
#   - thread / concurrency derivation        (concurrency.sh)
#   - per-sample resume / skip / force       (checkpoint.sh)
#   - real-time progress to web layer        (runtime.sh)
#   - graceful Ctrl+C with process-group kill (trap_add)
#   - structured event emission              (runtime.sh)
#   - multi-sample loop with FIFO semaphore  (concurrency.sh)
#   - per-sample retry on transient failure
#
# Usage from a pipeline.sh:
#
#   #!/usr/bin/env bash
#   source pipelines/lib/orchestrator.sh
#
#   run_preAlignmentQC_for_sample() {
#     local sample="$1"
#     # ... do work for one sample, emit progress, write log ...
#   }
#
#   pipeline_run preAlignmentQC
#
# Required environment / variables (typically loaded from ConfigFile):
#   work_dir            – derived workdir (created by ngs_build_workdir)
#   ntask_per_run       – concurrency hint (integer or "ALL")
#   samples[]           – populated by ngs_load_sample_info
#   force               – TRUE to re-run completed samples (default: FALSE)
#   retry               – number of retries per sample on failure (default: 1)
###############################################################################

# Locate the lib directory robustly so steps can source us from anywhere.
_NGS_LIB_DIR="${_NGS_LIB_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export _NGS_LIB_DIR

# shellcheck source=runtime.sh
source "$_NGS_LIB_DIR/runtime.sh"
# shellcheck source=concurrency.sh
source "$_NGS_LIB_DIR/concurrency.sh"
# shellcheck source=checkpoint.sh
source "$_NGS_LIB_DIR/checkpoint.sh"
# shellcheck source=workdir.sh
source "$_NGS_LIB_DIR/workdir.sh"
# shellcheck source=igenomes.sh
source "$_NGS_LIB_DIR/igenomes.sh"

# `trap_add` may be provided by the legacy GlobalFunction.sh; if not,
# fall back to a single SIGINT/SIGTERM trap that kills the process group.
if ! declare -f trap_add >/dev/null 2>&1; then
  trap_add() {
    local cmd="$1"; shift
    trap "$cmd" "$@"
  }
fi

# ---------------------------------------------------------------------------
# pipeline_run <pipeline_name>
#
# Iterates over $samples, dispatching `run_<pipeline_name>_for_sample` for
# each sample under FIFO-controlled concurrency.
#
# A sample is skipped when:
#   - $force != "TRUE", and
#   - $work_dir/$sample/.ngs_complete.<pipeline> exists
#
# Failed samples are retried up to $retry times; persistent failures are
# recorded to $work_dir/.ngs_failed.<pipeline> for the web UI to surface.
# ---------------------------------------------------------------------------
pipeline_run() {
  local pipeline="${1:?pipeline_run: name required}"
  local fn="run_${pipeline}_for_sample"
  if ! declare -f "$fn" >/dev/null 2>&1; then
    log_error "pipeline_run: function $fn is not defined"
    return 1
  fi

  : "${work_dir:?work_dir is not set}"
  : "${ntask_per_run:=1}"
  : "${force:=FALSE}"
  : "${retry:=1}"
  if (( ${#samples[@]} == 0 )); then
    log_error "pipeline_run: samples[] is empty (did you call ngs_load_sample_info?)"
    return 1
  fi

  ngs_compute_threads "${#samples[@]}"
  log_info "pipeline=$pipeline samples=${#samples[@]} concurrency=$ntask_per_run threads=$threads"
  emit_status "Starting $pipeline on ${#samples[@]} samples"
  emit_progress 0

  ngs_fifo_open "$ntask_per_run"

  local failed_log="$work_dir/.ngs_failed.$pipeline"
  : > "$failed_log"

  local total="${#samples[@]}"
  local done_count=0
  local pids=()

  local sample
  for sample in "${samples[@]}"; do
    ngs_fifo_acquire
    {
      local sample_dir="$work_dir/$sample"
      mkdir -p "$sample_dir"

      if [[ "$force" == "TRUE" ]]; then
        ngs_pipeline_clear_complete "$sample_dir" "$pipeline"
      fi

      if ngs_pipeline_completed "$sample_dir" "$pipeline"; then
        log_info "[$sample] skip (already complete)"
        emit_status "[$sample] skip"
      else
        local attempt=1 rc=1
        while (( attempt <= retry )) && (( rc != 0 )); do
          if (( attempt > 1 )); then
            log_warn "[$sample] retry attempt $attempt/$retry"
            emit_status "[$sample] retry $attempt/$retry"
          fi
          # Run the sample's pipeline body. Tolerate failure so we can
          # log + continue with other samples.
          set +e
          "$fn" "$sample"
          rc=$?
          set -e
          attempt=$((attempt + 1))
        done
        if (( rc == 0 )); then
          ngs_pipeline_mark_complete "$sample_dir" "$pipeline"
          log_ok "[$sample] $pipeline complete"
        else
          log_error "[$sample] $pipeline FAILED after $retry attempt(s)"
          printf '%s\n' "$sample" >> "$failed_log"
        fi
      fi

      # Best-effort progress: percentage of completed/skipped samples.
      # We use a flock'd counter file so concurrent samples don't race.
      (
        flock 9
        done_count=$(($(cat "$work_dir/.ngs_done.$pipeline" 2>/dev/null || echo 0) + 1))
        printf '%s' "$done_count" > "$work_dir/.ngs_done.$pipeline"
      ) 9>"$work_dir/.ngs_done.lock.$pipeline"
      done_count="$(cat "$work_dir/.ngs_done.$pipeline" 2>/dev/null || echo 0)"
      emit_progress "$(( 100 * done_count / total ))"

      ngs_fifo_release
    } &
    pids+=("$!")
  done

  # Wait for everything; allow Ctrl+C to propagate via the process-group
  # trap installed by ngs_fifo_open.
  wait

  ngs_fifo_close
  rm -f "$work_dir/.ngs_done.$pipeline" "$work_dir/.ngs_done.lock.$pipeline"

  local failed=0
  if [[ -s "$failed_log" ]]; then
    failed=$(wc -l < "$failed_log")
  fi
  emit_progress 100
  if (( failed > 0 )); then
    emit_status "$pipeline finished with $failed failed sample(s)"
    log_error "Pipeline $pipeline: $failed failed sample(s) — see $failed_log"
    return 1
  fi
  emit_status "$pipeline finished cleanly"
  log_ok "Pipeline $pipeline: ${#samples[@]} samples succeeded"
  return 0
}

# ---------------------------------------------------------------------------
# Helper: read a YAML metadata file (very small subset, top-level scalars
# and `default_params:` map). Pipelines call this to load their meta.yml
# without an external YAML parser.
#
# Usage:
#   declare -A meta
#   pipeline_meta_load /path/to/meta.yml meta
#   echo "${meta[name]}"          # → preAlignmentQC
#   echo "${meta[default_params.min_quality]}"  # → 20
# ---------------------------------------------------------------------------
pipeline_meta_load() {
  local file="${1:?meta yml path required}"
  local -n _out="$2"
  local in_params=0 line key val
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    if [[ "$line" =~ ^default_params: ]]; then
      in_params=1
      continue
    fi
    if (( in_params )) && [[ "$line" =~ ^[[:space:]]+([A-Za-z_][A-Za-z_0-9]*):[[:space:]]*(.*)$ ]]; then
      key="${BASH_REMATCH[1]}"
      val="${BASH_REMATCH[2]}"
      val="${val%\"}"; val="${val#\"}"; val="${val%\'}"; val="${val#\'}"
      _out["default_params.$key"]="$val"
      continue
    fi
    # leaving the params block when indentation drops
    if (( in_params )) && [[ "$line" =~ ^[A-Za-z] ]]; then
      in_params=0
    fi
    if [[ "$line" =~ ^([A-Za-z_][A-Za-z_0-9]*):[[:space:]]*(.*)$ ]]; then
      key="${BASH_REMATCH[1]}"
      val="${BASH_REMATCH[2]}"
      val="${val%\"}"; val="${val#\"}"; val="${val%\'}"; val="${val#\'}"
      [[ -n "$val" ]] && _out["$key"]="$val"
    fi
  done < "$file"
}
