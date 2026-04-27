#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/profiles.sh
#
# Execution profiles. A profile decides *where* a per-sample stage runs:
#   - local   → current shell (default; matches existing behaviour)
#   - slurm   → sbatch --wait wrapper
#   - sge     → qsub -sync y wrapper
#   - lsf     → bsub -K wrapper
#
# Selected via NGS_PROFILE env (or `ngsmodule run --profile slurm`).
# Resources come from the pipeline's meta.yml `resources:` block:
#
#   resources:
#     mem: 32G
#     time: 4h
#     cpus_per_task: 8
#     queue: high-mem      # SGE/LSF queue name
#     partition: short     # Slurm partition
#     extra: "--gres=gpu:1"# additional native flags
#
# Each profile script must define `profile_submit_sample <pipeline> <sample>`.
# pipeline_run_one() calls it instead of running the per-sample function
# directly when NGS_PROFILE != "local".
###############################################################################

# Read the resources block from meta.yml into an associative array.
# Output keys: mem, time, cpus_per_task, queue, partition, extra.
pipeline_resources() {
  local pipeline="${1:?pipeline required}"
  local -n _out="${2:?out array required}"
  local meta="${PIPELINES_DIR:-$NGSMODULE_ROOT/pipelines/core}/$pipeline/meta.yml"
  [[ ! -f "$meta" ]] && return 0

  local in_block=0 line key val
  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%$'\r'}"
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    if [[ "$line" =~ ^resources:[[:space:]]*$ ]]; then
      in_block=1; continue
    fi
    if (( in_block )) && [[ "$line" =~ ^[A-Za-z_] ]]; then
      in_block=0
    fi
    (( in_block )) || continue
    if [[ "$line" =~ ^[[:space:]]+([A-Za-z_][A-Za-z_0-9]*):[[:space:]]*(.*)$ ]]; then
      key="${BASH_REMATCH[1]}"
      val="${BASH_REMATCH[2]}"
      val="${val%\"}"; val="${val#\"}"
      val="${val%\'}"; val="${val#\'}"
      _out["$key"]="$val"
    fi
  done < "$meta"
}

# Load the configured profile (defaulting to local).
# Idempotent.
profile_load() {
  local name="${NGS_PROFILE:-local}"
  local file="${NGSMODULE_ROOT:-/}/pipelines/profiles/${name}.sh"
  if [[ ! -f "$file" ]]; then
    echo "ERROR: unknown profile '$name' (no such file: $file)" >&2
    return 1
  fi
  if ! declare -f profile_submit_sample >/dev/null 2>&1; then
    # shellcheck source=/dev/null
    source "$file"
  fi
  return 0
}

# Wrapper: call profile_submit_sample if a profile is loaded; else run inline.
# Pipelines should not call this directly — pipeline_run_one already does.
profile_dispatch() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"
  local fn="run_${pipeline}_for_sample"
  if [[ "${NGS_PROFILE:-local}" == "local" ]]; then
    "$fn" "$sample"
    return $?
  fi
  profile_load || return $?
  profile_submit_sample "$pipeline" "$sample"
}
