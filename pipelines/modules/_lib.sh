#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/_lib.sh
#
# Shared module helpers. Pipelines source individual tool modules via
# `module_load <tool>`; modules source this file to gain access to:
#
#   module_load <tool>          # source pipelines/modules/<tool>.sh once
#   module_run "<desc>" cmd...  # run_step + auto-track tool usage
#   module_tools_json           # emit the JSON array of tools used so far
#                               # in this pipeline run (for state.json)
#
# Inspired by nf-core/modules but: pure bash, no DSL, no per-process
# isolation. The contract is just: "a function that takes named args
# and writes to a known location". Pipelines compose modules by calling
# their entry-point functions in sequence.
###############################################################################

# Already loaded? bail.
[[ -n "${_NGS_MODULES_LIB_LOADED:-}" ]] && return 0
_NGS_MODULES_LIB_LOADED=1

_MODULES_DIR="${_MODULES_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export _MODULES_DIR

# Tools observed this session (fed to ngs_state_tools_json by the pipeline).
declare -gA _MODULES_TOOLS_USED

# ---------------------------------------------------------------------------
# Source a module by name. Idempotent — safe to call multiple times.
# Errors out clearly if the module file is missing.
# ---------------------------------------------------------------------------
module_load() {
  local tool="${1:?module_load: tool name required}"
  local file="$_MODULES_DIR/${tool}.sh"
  if [[ ! -f "$file" ]]; then
    echo "ERROR: module not found: $file" >&2
    return 1
  fi
  # Guard against double-sourcing.
  local guard="_NGS_MOD_${tool//-/_}_LOADED"
  if [[ -z "${!guard:-}" ]]; then
    # shellcheck source=/dev/null
    source "$file"
    printf -v "$guard" '%s' 1
    export "${guard?}"
  fi
}

# ---------------------------------------------------------------------------
# Thin run_step wrapper that also marks the tool as used so we can build
# a state.json `tools` array later. Pass the tool binary as first cmd arg
# so we can sniff its name; or set the env var NGS_MODULE_CURRENT_TOOL
# explicitly before calling.
# ---------------------------------------------------------------------------
module_run() {
  local desc="${1:?module_run: description required}"; shift
  local tool="${NGS_MODULE_CURRENT_TOOL:-${1##*/}}"
  _MODULES_TOOLS_USED["$tool"]=1
  if declare -f is_dry_run >/dev/null && is_dry_run; then
    printf '::dry-run::%s\n' "$*"
    return 0
  fi
  if declare -f run_step >/dev/null; then
    run_step "$desc" "$@"
  else
    "$@"
  fi
}

# ---------------------------------------------------------------------------
# Emit a JSON array of all tools used by this pipeline run. Pipelines hand
# this to ngs_state_stage_end --tools.
# ---------------------------------------------------------------------------
module_tools_json() {
  local args=()
  local t
  for t in "${!_MODULES_TOOLS_USED[@]}"; do
    args+=("$t")
  done
  if declare -f ngs_state_tools_json >/dev/null; then
    ngs_state_tools_json "${args[@]}"
  else
    printf '['
    local first=1
    for t in "${args[@]}"; do
      (( first )) || printf ','
      printf '{"name":"%s"}' "$t"
      first=0
    done
    printf ']'
  fi
}

# Reset between pipelines (when one shell runs multiple).
module_reset() { _MODULES_TOOLS_USED=(); }
