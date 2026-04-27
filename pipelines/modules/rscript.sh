#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/modules/rscript.sh — Run R scripts with package presence check.
#
# Functions:
#   module_rscript_check_packages <pkg1> <pkg2>...
#       Verifies that each named R package is installed. Returns 1 with
#       a clear error if any are missing. Honours $NGS_RSCRIPT_BIN.
#
#   module_rscript_run --script <path.R> --args "<arg1 arg2 ...>"
#                      [--cwd <dir>] [--required <pkg1,pkg2,...>]
#       Optionally pre-checks required packages, then runs the script
#       under module_run so the call is benchmarked + tool-tracked.
#
# Honoured env vars:
#   NGS_RSCRIPT_BIN — path to Rscript binary (default: `Rscript`)
###############################################################################

source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

_rscript_bin() { echo "${NGS_RSCRIPT_BIN:-Rscript}"; }

module_rscript_check_packages() {
  if (( $# == 0 )); then return 0; fi
  local rs
  rs="$(_rscript_bin)"
  if ! command -v "$rs" >/dev/null 2>&1; then
    echo "ERROR: $rs not on PATH (set NGS_RSCRIPT_BIN to override)" >&2
    return 1
  fi
  local missing
  # Single Rscript invocation that prints space-separated missing names.
  # We intentionally use double quotes inside the R one-liner so $@ expands.
  missing="$("$rs" --vanilla -e "
    pkgs <- c($(printf '\"%s\",' "$@" | sed 's/,$//'))
    miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (length(miss) > 0) cat(paste(miss, collapse = ' '))
  " 2>/dev/null)"
  if [[ -n "$missing" ]]; then
    echo "ERROR: missing R package(s): $missing" >&2
    echo "       install via: BiocManager::install(c($(echo "$missing" | sed "s/[^ ]*/\'&\'/g; s/ /,/g")))" >&2
    return 1
  fi
  return 0
}

module_rscript_run() {
  local script="" cwd="" required=""
  local -a script_args=()
  while (( $# > 0 )); do
    case "$1" in
      --script)   script="$2"; shift 2 ;;
      --cwd)      cwd="$2"; shift 2 ;;
      --required) required="$2"; shift 2 ;;
      --args)     read -ra script_args <<< "$2"; shift 2 ;;
      *) echo "module_rscript_run: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${script:?--script required}"
  if [[ ! -f "$script" ]] && ! { declare -f is_dry_run >/dev/null && is_dry_run; }; then
    echo "ERROR: Rscript path not found: $script" >&2
    return 1
  fi

  if [[ -n "$required" ]]; then
    # Comma- or space-separated list of required packages.
    local IFS=', '
    # shellcheck disable=SC2206
    local pkgs=($required)
    if ! module_rscript_check_packages "${pkgs[@]}"; then return 1; fi
  fi

  local rs
  rs="$(_rscript_bin)"

  if [[ -n "$cwd" ]]; then
    NGS_MODULE_CURRENT_TOOL=Rscript \
      module_run "Rscript $(basename "$script")" \
        bash -c "cd '$cwd' && '$rs' --vanilla '$script' $(printf "%q " "${script_args[@]}")"
  else
    NGS_MODULE_CURRENT_TOOL=Rscript \
      module_run "Rscript $(basename "$script")" \
        "$rs" --vanilla "$script" "${script_args[@]}"
  fi
}
