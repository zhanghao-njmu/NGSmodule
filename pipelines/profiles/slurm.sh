#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/profiles/slurm.sh — submit each sample × stage as one Slurm job.
#
# We use `sbatch --wait` so the orchestrator's per-sample loop blocks
# (and the FIFO semaphore still throttles concurrent submissions, which
# is sometimes what you want; otherwise set ntask=ALL and let Slurm
# manage the queue).
#
# Job script convention:
#   $work_dir/$sample/.slurm/${pipeline}.sh
#   $work_dir/$sample/.slurm/${pipeline}.out
#   $work_dir/$sample/.slurm/${pipeline}.err
#
# Resources are read from meta.yml — see pipelines/lib/profiles.sh.
###############################################################################

profile_submit_sample() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"

  : "${work_dir:?work_dir not set}"
  : "${NGSMODULE_ROOT:?NGSMODULE_ROOT not set}"

  local jdir="$work_dir/$sample/.slurm"
  mkdir -p "$jdir"
  local script="$jdir/${pipeline}.sh"
  local out_log="$jdir/${pipeline}.out"
  local err_log="$jdir/${pipeline}.err"

  declare -A R=()
  pipeline_resources "$pipeline" R

  # Render a self-contained job script. We re-source ConfigFile inside
  # so the job inherits a clean environment (matches Slurm's --export=NONE
  # convention without forcing it).
  cat >"$script" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=ngs.${pipeline}.${sample}
#SBATCH --output=${out_log}
#SBATCH --error=${err_log}
${R[partition]:+#SBATCH --partition=\"${R[partition]}\"}
${R[mem]:+#SBATCH --mem=\"${R[mem]}\"}
${R[time]:+#SBATCH --time=\"${R[time]}\"}
${R[cpus_per_task]:+#SBATCH --cpus-per-task=${R[cpus_per_task]}}
${R[extra]:+#SBATCH ${R[extra]}}
set -euo pipefail
cd "${NGSMODULE_ROOT}"
${ConfigFile:+source "${ConfigFile}"}
export NGS_PROFILE=local   # inside the node, run inline (no recursion)
export work_dir="${work_dir}"
${SampleInfoFile:+export SampleInfoFile="${SampleInfoFile}"}
source pipelines/core/${pipeline}/pipeline.sh
${SampleInfoFile:+ngs_load_sample_info "${SampleInfoFile}"}
run_${pipeline}_for_sample "${sample}"
EOF
  chmod +x "$script"

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not on PATH; cannot run profile=slurm" >&2
    return 1
  fi

  log_info "[$sample] $pipeline → sbatch $script"
  sbatch --wait "$script" >>"$out_log" 2>>"$err_log"
}
