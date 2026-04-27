#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/runtime.sh
#
# Shared runtime helpers for NGSmodule pipeline scripts. Source this file at
# the top of any GeneralSteps/Analysis script to get:
#
#   - `set -euo pipefail` with descriptive failure traps
#   - Structured progress events (`::progress::N`) consumed by the web worker
#   - `::artifact::PATH` and `::metric::KEY=VALUE` events for indexing results
#   - Color-coded log levels (`log_info`, `log_warn`, `log_error`, `log_ok`)
#   - `--dry-run` short-circuit (set NGS_DRY_RUN=1 or pass `--dry-run`)
#   - `run_step "label" cmd args...` wrapper that emits begin/end events and
#     respects dry-run mode
#
# Backwards-compatible: scripts that don't source this file behave exactly as
# before. Existing scripts can opt in incrementally.
###############################################################################

# ---------------------------------------------------------------------------
# Strict mode (only when explicitly opted in via NGS_STRICT=1).
# We don't enforce this on every script yet because legacy postcheck blocks
# rely on continuing after individual tool failures.
# ---------------------------------------------------------------------------
if [[ "${NGS_STRICT:-0}" == "1" ]]; then
  set -euo pipefail
  # On unexpected exit, surface the failing line + command to stderr so the
  # web worker captures it in the task log.
  trap 'rc=$?; echo "::error::line=${LINENO} cmd=${BASH_COMMAND} rc=${rc}" >&2; exit $rc' ERR
fi

# ---------------------------------------------------------------------------
# Dry-run support.
#
# When NGS_DRY_RUN=1 the wrapper prints commands instead of executing them.
# Useful for the web UI to preview a pipeline before committing CPU time.
# ---------------------------------------------------------------------------
NGS_DRY_RUN="${NGS_DRY_RUN:-0}"

is_dry_run() { [[ "$NGS_DRY_RUN" == "1" ]]; }

# ---------------------------------------------------------------------------
# Structured event protocol.
#
# These markers are emitted on stdout. The Celery worker tails the log file
# and parses them. They are also human-readable in the captured log.
# ---------------------------------------------------------------------------

# Emit a progress percentage (0-100). Workers convert this to a
# float in [0.0, 1.0] for the realtime WebSocket.
emit_progress() {
  local pct="$1"
  printf '::progress::%s\n' "$pct"
}

# Emit a status message visible in the UI's task pane.
emit_status() {
  local msg="$*"
  printf '::status::%s\n' "$msg"
}

# Register a result artefact (file path) so the web layer can list it
# under Results without filesystem scanning.
emit_artifact() {
  local kind="$1"; shift
  local path="$1"; shift
  printf '::artifact::kind=%s path=%s\n' "$kind" "$path"
}

# Register a metric (e.g. "total_reads=12345678" or "duplication_rate=0.18").
emit_metric() {
  local key="$1"; shift
  local value="$1"; shift
  printf '::metric::%s=%s\n' "$key" "$value"
}

# ---------------------------------------------------------------------------
# Colored logging (only when stdout is a TTY; web workers see clean text).
# ---------------------------------------------------------------------------
if [[ -t 1 ]]; then
  _C_RESET=$'\033[0m'; _C_RED=$'\033[31m'; _C_GREEN=$'\033[32m'
  _C_YELLOW=$'\033[33m'; _C_BLUE=$'\033[34m'; _C_BOLD=$'\033[1m'
else
  _C_RESET=''; _C_RED=''; _C_GREEN=''
  _C_YELLOW=''; _C_BLUE=''; _C_BOLD=''
fi

_log() {
  local prefix="$1"; shift
  printf '%s[%s]%s %s\n' "$prefix" "$(date '+%Y-%m-%d %H:%M:%S')" "$_C_RESET" "$*"
}

log_info()  { _log "$_C_BLUE${_C_BOLD}INFO" "$@"; }
log_ok()    { _log "$_C_GREEN${_C_BOLD}OK"   "$@"; }
log_warn()  { _log "$_C_YELLOW${_C_BOLD}WARN" "$@" >&2; }
log_error() { _log "$_C_RED${_C_BOLD}ERROR"  "$@" >&2; }

# ---------------------------------------------------------------------------
# `run_step "label" cmd args...`
#
# Wraps a tool invocation with begin/end events and respects dry-run.
# Returns the underlying exit code so callers can branch.
#
# Example:
#   run_step "fastp trim" fastp -i s.fq.gz -o s.trim.fq.gz --thread 8
# ---------------------------------------------------------------------------
run_step() {
  local label="$1"; shift
  log_info "→ ${label}"
  emit_status "${label}"
  if is_dry_run; then
    printf '::dry-run::%s\n' "$*"
    return 0
  fi
  "$@"
  local rc=$?
  if [[ $rc -eq 0 ]]; then
    log_ok "✓ ${label}"
  else
    log_error "✗ ${label} (rc=$rc)"
  fi
  return $rc
}

# ---------------------------------------------------------------------------
# Resolve a task id from environment (set by the web worker via
# `env={"NGS_TASK_ID": ...}`). Empty string when run from CLI.
# ---------------------------------------------------------------------------
ngs_task_id() { printf '%s' "${NGS_TASK_ID:-}"; }

# ---------------------------------------------------------------------------
# Quick required-tool check. Outputs each tool's resolved path + version (if
# obtainable) at the start of a run so reproducibility is captured in the log.
#
# Usage: ngs_require fastp samtools STAR
# ---------------------------------------------------------------------------
ngs_require() {
  local missing=0
  for tool in "$@"; do
    if command -v "$tool" >/dev/null 2>&1; then
      local ver
      ver="$("$tool" --version 2>&1 | head -1 || echo unknown)"
      log_info "tool: $tool ($(command -v "$tool"))  $ver"
    else
      log_error "tool not found in PATH: $tool"
      missing=$((missing + 1))
    fi
  done
  if [[ $missing -gt 0 ]]; then
    log_error "$missing required tools missing — aborting"
    exit 127
  fi
}
