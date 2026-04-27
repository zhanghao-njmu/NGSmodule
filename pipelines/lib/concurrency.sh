#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/concurrency.sh
#
# Per-sample concurrency control via a kernel FIFO semaphore.
#
# Why FIFO and not GNU parallel / xargs?
#   - Zero external dependencies (any POSIX bash).
#   - Survives Ctrl+C cleanly via process-group kill (`kill -- -$$`).
#   - Composable with `wait` so the parent shell knows when everything
#     finished without polling.
#
# Public API:
#   ngs_compute_threads <total_tasks>   # populates threads, threads_fastp,
#                                       # threads_featurecounts, threads_bismark
#                                       # based on ntask_per_run + total_threads
#                                       # (auto-detected from /proc/cpuinfo)
#   ngs_fifo_open <ntask_per_run>       # opens FD 1000 with N tokens
#   ngs_fifo_acquire                    # blocks until a token is available
#   ngs_fifo_release                    # returns a token to the pool
#   ngs_fifo_close                      # cleans up FIFO + FD
#
# The classic per-sample loop becomes:
#
#   ngs_compute_threads "${#samples[@]}"
#   ngs_fifo_open "$ntask_per_run"
#   for sample in "${samples[@]}"; do
#     ngs_fifo_acquire
#     {
#       run_pipeline_for_sample "$sample"
#       ngs_fifo_release
#     } &
#   done
#   wait
#   ngs_fifo_close
#
# Variables that get exported so step scripts can use them directly:
#   ntask_per_run, total_threads, threads, threads_fastp,
#   threads_featurecounts, threads_bismark
###############################################################################

# ---------------------------------------------------------------------------
# Total host threads. Use NPROC if set (lets caller restrict CPUs in cgroup
# environments); otherwise count all CPUs from /proc.
# ---------------------------------------------------------------------------
ngs_total_threads() {
  if [[ -n "${NPROC:-}" ]]; then
    printf '%s' "$NPROC"
    return
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc
  else
    grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 4
  fi
}

# ---------------------------------------------------------------------------
# Resolve `ntask_per_run` (user-set; can be integer or "ALL") into a number
# and derive per-task thread counts. The thread caps mirror what was in
# LoadConfig.sh — they reflect real-world tooling sweet spots:
#
#   - 32 threads is the upper limit before NUMA cross-domain becomes painful
#     for STAR / BWA-MEM2.
#   - fastp is I/O-bound past 16 threads.
#   - featureCounts scales to 64 threads on big BAMs.
#   - bismark spawns 8 sub-threads per mapper, so divide by 8.
#
# Args: $1 = total_tasks (number of samples this run)
# Globals consumed: ntask_per_run
# Globals produced (exported): total_threads, threads, threads_fastp,
#                              threads_featurecounts, threads_bismark
# ---------------------------------------------------------------------------
ngs_compute_threads() {
  local total_task="${1:-1}"
  total_threads="$(ngs_total_threads)"
  export total_threads

  if [[ -z "${ntask_per_run:-}" ]]; then
    ntask_per_run=1
  fi

  if [[ "$ntask_per_run" =~ ^[0-9]+$ ]]; then
    : # already numeric
  elif [[ "$ntask_per_run" == "ALL" ]]; then
    if (( total_task > total_threads )); then
      ntask_per_run=$total_threads
    else
      ntask_per_run=$total_task
    fi
  else
    echo "ERROR: ntask_per_run must be 'ALL' or an integer (got: $ntask_per_run)" >&2
    return 1
  fi
  export ntask_per_run

  threads=$(( total_threads / ntask_per_run ))
  (( threads < 1 )) && threads=1
  (( threads > 32 )) && threads=32
  export threads

  if (( threads > 16 )); then
    threads_fastp=16
  else
    threads_fastp=$threads
  fi
  export threads_fastp

  if (( threads > 64 )); then
    threads_featurecounts=64
  else
    threads_featurecounts=$threads
  fi
  export threads_featurecounts

  threads_bismark=$(( threads / 8 ))
  (( threads_bismark < 1 )) && threads_bismark=1
  export threads_bismark
}

# ---------------------------------------------------------------------------
# FIFO-based semaphore. We use FD 1000 to avoid colliding with anything
# user scripts might open. The mkfifo + immediate rm trick keeps a zombie
# file out of the cwd while the FD stays open.
# ---------------------------------------------------------------------------
ngs_fifo_open() {
  local n="${1:?ngs_fifo_open: ntask_per_run required}"
  local fifo
  fifo="${TMPDIR:-/tmp}/ngs_fifo_$$.$RANDOM"
  mkfifo "$fifo"
  exec 1000<>"$fifo"
  rm -f "$fifo"
  local i
  for (( i = 1; i <= n; i++ )); do
    echo >&1000
  done
  # Best-effort cleanup if the script is killed.
  if declare -f trap_add >/dev/null 2>&1; then
    trap_add 'ngs_fifo_close' SIGINT SIGTERM EXIT
  else
    trap 'ngs_fifo_close' SIGINT SIGTERM EXIT
  fi
}

ngs_fifo_acquire() {
  local _token
  read -r -u 1000 _token
}

ngs_fifo_release() {
  echo >&1000
}

ngs_fifo_close() {
  exec 1000>&- 2>/dev/null || true
  exec 1000<&- 2>/dev/null || true
}
