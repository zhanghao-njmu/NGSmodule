#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/profiles/lsf.sh — IBM Spectrum LSF submission.
#
# Uses `bsub -K` so the orchestrator blocks per submission.
# Resources translated:
#   mem            → -M (per-task memory limit, MB)
#   time           → -W (HH:MM)
#   cpus_per_task  → -n
#   queue          → -q
#   extra          → appended verbatim
###############################################################################

profile_submit_sample() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"
  : "${work_dir:?work_dir not set}"

  local jdir="$work_dir/$sample/.lsf"
  mkdir -p "$jdir"
  local script="$jdir/${pipeline}.sh"
  local out_log="$jdir/${pipeline}.out"
  local err_log="$jdir/${pipeline}.err"

  declare -A R=()
  pipeline_resources "$pipeline" R

  cat >"$script" <<EOF
#!/usr/bin/env bash
#BSUB -J ngs.${pipeline}.${sample}
#BSUB -o ${out_log}
#BSUB -e ${err_log}
${R[queue]:+#BSUB -q "${R[queue]}"}
${R[mem]:+#BSUB -M "${R[mem]}"}
${R[time]:+#BSUB -W "${R[time]}"}
${R[cpus_per_task]:+#BSUB -n ${R[cpus_per_task]}}
${R[extra]:+#BSUB ${R[extra]}}
set -euo pipefail
cd "${NGSMODULE_ROOT}"
${ConfigFile:+source "${ConfigFile}"}
export NGS_PROFILE=local
export work_dir="${work_dir}"
${SampleInfoFile:+export SampleInfoFile="${SampleInfoFile}"}
source pipelines/core/${pipeline}/pipeline.sh
${SampleInfoFile:+ngs_load_sample_info "${SampleInfoFile}"}
run_${pipeline}_for_sample "${sample}"
EOF
  chmod +x "$script"

  if ! command -v bsub >/dev/null 2>&1; then
    echo "ERROR: bsub not on PATH; cannot run profile=lsf" >&2
    return 1
  fi
  log_info "[$sample] $pipeline → bsub $script"
  bsub -K <"$script" >>"$out_log" 2>>"$err_log"
}
