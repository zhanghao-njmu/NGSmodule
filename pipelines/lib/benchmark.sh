#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/benchmark.sh
#
# Capture wall-clock + max RSS + CPU% per command using GNU /usr/bin/time -v,
# then aggregate into a per-stage benchmark blob that lands in state.json.
#
# Inspired by Snakemake's `benchmark:` directive but framed around our
# sample × stage state model. Aggregation is "max" for resource peaks and
# "sum" for wall time so the final number reflects the whole stage.
#
# Usage from the orchestrator:
#
#   ngs_benchmark_begin "$sample" "$pipeline"
#   ... (run_step calls automatically wrap with /usr/bin/time -v) ...
#   ngs_benchmark_json   # prints {"max_rss_kb":N,"cpu_pct":P,"wall_seconds":S}
#   ngs_benchmark_end
#
# Disabled when /usr/bin/time is missing or NGS_NO_BENCH=1.
###############################################################################

# Runtime state — global so run_step in any process inherits the file path
# via env. We intentionally use a temp file (one append per command) and
# aggregate only at end, so concurrent run_step invocations within one
# stage don't race on a counter.
declare -g _NGS_BENCH_FILE=""
declare -g _NGS_BENCH_TIME_BIN=""

ngs_benchmark_available() {
  [[ "${NGS_NO_BENCH:-0}" != "1" ]] && [[ -x /usr/bin/time ]]
}

# Begin per-stage capture. No-op if disabled.
ngs_benchmark_begin() {
  ngs_benchmark_available || { _NGS_BENCH_FILE=""; return 0; }
  local sample="${1:?sample required}"
  local pipeline="${2:?pipeline required}"
  _NGS_BENCH_TIME_BIN=/usr/bin/time
  _NGS_BENCH_FILE="$(mktemp -t ngs-bench.XXXXXX)"
  export _NGS_BENCH_FILE
}

# Wrap a command with /usr/bin/time -v, appending the parsed report to
# the per-stage file. Called by run_step when _NGS_BENCH_FILE is set.
ngs_benchmark_wrap() {
  [[ -z "${_NGS_BENCH_FILE:-}" ]] && { "$@"; return $?; }
  local tmp
  tmp="$(mktemp -t ngs-bench-step.XXXXXX)"
  "$_NGS_BENCH_TIME_BIN" -v -o "$tmp" "$@"
  local rc=$?
  # Append normalised key=value lines to the stage file.
  if [[ -s "$tmp" ]]; then
    awk -F': ' '
      /Maximum resident set size/ {printf "max_rss_kb=%s\n", $2}
      /Percent of CPU/            {gsub("%","",$2); printf "cpu_pct=%s\n", $2}
      /Elapsed \(wall clock\)/    {printf "wall_clock=%s\n", $2}
      /User time/                 {printf "user_seconds=%s\n", $2}
      /System time/               {printf "sys_seconds=%s\n", $2}
    ' "$tmp" >> "$_NGS_BENCH_FILE"
    echo "---" >> "$_NGS_BENCH_FILE"
  fi
  rm -f "$tmp"
  return $rc
}

# Convert "h:mm:ss" / "m:ss[.frac]" to seconds. Pure awk for portability.
_ngs_bench_secs() {
  awk -v t="$1" 'BEGIN {
    n = split(t, a, ":")
    if (n == 3) print a[1]*3600 + a[2]*60 + a[3]
    else if (n == 2) print a[1]*60 + a[2]
    else print t + 0
  }'
}

# Emit a JSON object for the stage. Empty object when capture disabled or
# nothing was run.
ngs_benchmark_json() {
  [[ -z "${_NGS_BENCH_FILE:-}" || ! -s "$_NGS_BENCH_FILE" ]] && {
    printf '{}'
    return 0
  }
  local max_rss=0 max_cpu=0 wall=0 user=0 sys=0
  local key val
  while IFS='=' read -r key val; do
    case "$key" in
      max_rss_kb)   (( val > max_rss )) && max_rss=$val ;;
      cpu_pct)
        # cpu_pct can be float — compare via awk, but our running max is int
        local v_int="${val%%.*}"
        (( v_int > max_cpu )) && max_cpu=$v_int
        ;;
      wall_clock)
        local s
        s=$(_ngs_bench_secs "$val")
        wall=$(awk -v a="$wall" -v b="$s" 'BEGIN{print a+b}')
        ;;
      user_seconds) user=$(awk -v a="$user" -v b="$val" 'BEGIN{print a+b}') ;;
      sys_seconds)  sys=$(awk -v a="$sys" -v b="$val" 'BEGIN{print a+b}') ;;
    esac
  done < "$_NGS_BENCH_FILE"
  printf '{"max_rss_kb":%s,"cpu_pct":%s,"wall_seconds":%s,"user_seconds":%s,"sys_seconds":%s}' \
    "$max_rss" "$max_cpu" "$wall" "$user" "$sys"
}

# End per-stage capture and clean up.
ngs_benchmark_end() {
  [[ -n "${_NGS_BENCH_FILE:-}" && -f "$_NGS_BENCH_FILE" ]] && rm -f "$_NGS_BENCH_FILE"
  _NGS_BENCH_FILE=""
}
