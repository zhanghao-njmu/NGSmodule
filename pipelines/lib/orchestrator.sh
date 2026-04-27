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
# shellcheck source=sample_state.sh
source "$_NGS_LIB_DIR/sample_state.sh"

# Pipelines registry — populated by the CLI dispatcher; orchestrator reads
# `requires:` from each pipeline's meta.yml to auto-resolve prerequisites.
NGSMODULE_ROOT="${NGSMODULE_ROOT:-$(cd "$_NGS_LIB_DIR/../.." && pwd)}"
export NGSMODULE_ROOT
PIPELINES_DIR="${PIPELINES_DIR:-$NGSMODULE_ROOT/pipelines/core}"

# ---------------------------------------------------------------------------
# Read `requires:` (a YAML list of prerequisite pipelines) from meta.yml.
# Returns the list space-separated on stdout. Empty if no requires block.
#
#   requires:
#     - preAlignmentQC
#     - postAlignmentQC
# ---------------------------------------------------------------------------
pipeline_requires() {
  local pipeline="$1"
  local meta="$PIPELINES_DIR/$pipeline/meta.yml"
  [[ ! -f "$meta" ]] && return 0
  awk '
    /^requires:/        {in_block=1; next}
    /^[A-Za-z_]+:/      {in_block=0}
    in_block && /^[[:space:]]*-/ {
      sub(/^[[:space:]]*-[[:space:]]*/, "")
      sub(/[[:space:]]*$/, "")
      print
    }
  ' "$meta"
}

# Read `version:` from meta.yml (default 1.0.0).
pipeline_version() {
  local pipeline="$1"
  local meta="$PIPELINES_DIR/$pipeline/meta.yml"
  [[ ! -f "$meta" ]] && { echo "0.0.0"; return; }
  local v
  v="$(awk '/^version:/ {print $2; exit}' "$meta")"
  v="${v//\"/}"; v="${v//\'/}"
  echo "${v:-1.0.0}"
}

# ---------------------------------------------------------------------------
# Run prerequisite pipelines for a sample if not already complete.
# Recursive: a prereq's own prereqs are resolved transitively.
#
# Honoured by `pipeline_run` unless NGS_NO_DEPS=1 (or --no-deps).
# ---------------------------------------------------------------------------
_resolved_prereqs=""
pipeline_resolve_prereqs() {
  local pipeline="$1" sample="$2"
  local prereq
  for prereq in $(pipeline_requires "$pipeline"); do
    # Skip if already complete
    if ngs_state_is_complete "$sample" "$prereq"; then
      continue
    fi
    # Mark visited so we don't loop on cyclic graphs
    case " $_resolved_prereqs " in *" $prereq "*) continue ;; esac
    _resolved_prereqs+=" $prereq"

    log_info "[$sample] missing prerequisite: $prereq → resolving"
    emit_status "[$sample] resolving prerequisite: $prereq"

    local prereq_script="$PIPELINES_DIR/$prereq/pipeline.sh"
    if [[ ! -f "$prereq_script" ]]; then
      log_error "[$sample] missing prerequisite '$prereq' has no pipeline.sh"
      return 1
    fi

    # Resolve transitively first
    pipeline_resolve_prereqs "$prereq" "$sample" || return 1

    # Run the prereq's per-sample function (must already be loaded since
    # the pipeline.sh sources orchestrator.sh which sources prereq).
    # We source the prereq's pipeline.sh in the current shell so its
    # function becomes callable.
    if ! declare -f "run_${prereq}_for_sample" >/dev/null 2>&1; then
      # shellcheck source=/dev/null
      source "$prereq_script"
    fi
    pipeline_run_one "$prereq" "$sample" || return 1
  done
  return 0
}

# ---------------------------------------------------------------------------
# Run a single sample through a single pipeline stage with full state
# tracking. Used both by pipeline_run (per-sample loop) and by
# pipeline_resolve_prereqs (recursive prerequisite resolution).
# ---------------------------------------------------------------------------
pipeline_run_one() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"
  local fn="run_${pipeline}_for_sample"
  if ! declare -f "$fn" >/dev/null 2>&1; then
    log_error "[$sample] $fn not defined for pipeline '$pipeline'"
    return 1
  fi

  if ngs_state_is_complete "$sample" "$pipeline" && [[ "${force:-FALSE}" != "TRUE" ]]; then
    log_info "[$sample] $pipeline already complete (state.json) → skip"
    return 0
  fi

  local sample_dir="$work_dir/$sample"
  mkdir -p "$sample_dir"
  ngs_state_init "$sample"

  local version
  version="$(pipeline_version "$pipeline")"
  ngs_state_stage_begin "$sample" "$pipeline" "$version"

  local rc=0
  local attempt=1
  while (( attempt <= ${retry:-1} )); do
    if (( attempt > 1 )); then
      log_warn "[$sample] $pipeline retry attempt $attempt/${retry:-1}"
      emit_status "[$sample] retry $attempt"
    fi
    set +e
    "$fn" "$sample"
    rc=$?
    set -e
    [[ $rc -eq 0 ]] && break
    attempt=$((attempt + 1))
  done

  if (( rc == 0 )); then
    # Step scripts may have populated state-end metadata via
    # ngs_state_stage_end already; if not, emit a minimal completion marker.
    if [[ "$(ngs_state_stage_status "$sample" "$pipeline")" != "completed" ]]; then
      ngs_state_stage_end "$sample" "$pipeline" completed
    fi
    ngs_state_emit_receipt "$sample" "$pipeline"
    log_ok "[$sample] $pipeline complete"
  else
    ngs_state_stage_end "$sample" "$pipeline" failed --error "exit code $rc after $((attempt - 1)) attempt(s)"
    log_error "[$sample] $pipeline FAILED"
  fi
  return $rc
}

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

  # Echo dependency graph so users can see what will actually run.
  local prereqs
  prereqs="$(pipeline_requires "$pipeline" | tr '\n' ' ')"
  if [[ -n "${prereqs// /}" ]] && [[ "${NGS_NO_DEPS:-0}" != "1" ]]; then
    log_info "pipeline=$pipeline requires: $prereqs"
  fi

  ngs_compute_threads "${#samples[@]}"
  log_info "pipeline=$pipeline samples=${#samples[@]} concurrency=$ntask_per_run threads=$threads"
  emit_status "Starting $pipeline on ${#samples[@]} samples"
  emit_progress 0

  ngs_fifo_open "$ntask_per_run"

  local failed_log="$work_dir/.ngs_failed.$pipeline"
  : > "$failed_log"

  local total="${#samples[@]}"

  local sample
  for sample in "${samples[@]}"; do
    ngs_fifo_acquire
    {
      # Resolve prerequisites first (recursive, transitive). Skipped
      # when NGS_NO_DEPS=1.
      if [[ "${NGS_NO_DEPS:-0}" != "1" ]]; then
        if ! pipeline_resolve_prereqs "$pipeline" "$sample"; then
          printf '%s\n' "$sample" >> "$failed_log"
          ngs_fifo_release
          exit 0
        fi
      fi

      if ! pipeline_run_one "$pipeline" "$sample"; then
        printf '%s\n' "$sample" >> "$failed_log"
      fi

      # Per-pipeline progress aggregation (flock-guarded counter).
      (
        flock 9
        local done_count
        done_count=$(($(cat "$work_dir/.ngs_done.$pipeline" 2>/dev/null || echo 0) + 1))
        printf '%s' "$done_count" > "$work_dir/.ngs_done.$pipeline"
      ) 9>"$work_dir/.ngs_done.lock.$pipeline"
      local done_count
      done_count="$(cat "$work_dir/.ngs_done.$pipeline" 2>/dev/null || echo 0)"
      emit_progress "$(( 100 * done_count / total ))"

      ngs_fifo_release
    } &
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
