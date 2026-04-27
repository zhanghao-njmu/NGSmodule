#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/profiles/sge.sh — Sun/Univa Grid Engine submission.
#
# Uses `qsub -sync y` so the orchestrator blocks per submission.
# Resources translated:
#   mem            → -l h_vmem
#   time           → -l h_rt
#   cpus_per_task  → -pe smp <N>
#   queue          → -q
#   extra          → appended verbatim
###############################################################################

profile_submit_sample() {
  local pipeline="${1:?pipeline required}"
  local sample="${2:?sample required}"
  : "${work_dir:?work_dir not set}"

  local jdir="$work_dir/$sample/.sge"
  mkdir -p "$jdir"
  local script="$jdir/${pipeline}.sh"
  local out_log="$jdir/${pipeline}.out"
  local err_log="$jdir/${pipeline}.err"

  declare -A R=()
  pipeline_resources "$pipeline" R

  cat >"$script" <<EOF
#!/usr/bin/env bash
#$ -N ngs.${pipeline}.${sample}
#$ -o ${out_log}
#$ -e ${err_log}
#$ -cwd
${R[queue]:+#$ -q "${R[queue]}"}
${R[mem]:+#$ -l h_vmem="${R[mem]}"}
${R[time]:+#$ -l h_rt="${R[time]}"}
${R[cpus_per_task]:+#$ -pe smp ${R[cpus_per_task]}}
${R[extra]:+#$ ${R[extra]}}
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

  if ! command -v qsub >/dev/null 2>&1; then
    echo "ERROR: qsub not on PATH; cannot run profile=sge" >&2
    return 1
  fi
  log_info "[$sample] $pipeline → qsub $script"
  qsub -sync y "$script" >>"$out_log" 2>>"$err_log"
}
