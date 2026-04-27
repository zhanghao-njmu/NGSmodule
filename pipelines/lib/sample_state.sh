#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/sample_state.sh
#
# Sample-first state machine. The framework's defining concept.
#
# Rationale:
#   Wet-lab biologists ask "did sample Tumor_S01 finish QC?" — not "did
#   rule preAlignmentQC fire?". Existing workflow managers (Snakemake,
#   Nextflow) are DAG-first: they cache by rule × inputs hash, and surface
#   completion as "rules executed". NGSmodule is sample-first: every
#   sample owns a `.state.json` journal that explicitly records which
#   stages have completed, when, with which tools and parameters.
#
# This buys us:
#   - `ngsmodule status` becomes a sample × stage grid, not a DAG render
#   - Lineage: trace any output back to the inputs + tool version that
#     produced it without re-parsing logs
#   - Resume = "what's the next missing stage for this sample?" instead
#     of "is this Merkle hash the same as last time?"
#   - Cross-pipeline dependency resolution is trivial: "does sample X
#     have stage Y in its state?" — pure JSON lookup, no DAG walk
#
# State file format (JSON; one per sample, at $work_dir/$sample/.state.json):
#
#   {
#     "sample_id": "Tumor_S01",
#     "created_at": "2026-04-27T10:20:00Z",
#     "layout": "PE",
#     "group": "Treatment",
#     "stages": {
#       "preAlignmentQC": {
#         "status": "completed",          # pending | running | completed | failed
#         "started_at": "2026-04-27T10:20:01Z",
#         "completed_at": "2026-04-27T10:23:14Z",
#         "elapsed_seconds": 193,
#         "version": "1.0.0",
#         "params": {"min_quality": 20, "min_length": 50},
#         "tools": [
#           {"name": "fastp", "version": "0.23.4", "path": "/opt/conda/bin/fastp"}
#         ],
#         "inputs":  [{"path": "...", "sha256": "abc..."}],
#         "outputs": [{"path": "...", "sha256": "def...", "kind": "trimmed_fastq"}]
#       },
#       "Alignment": { "status": "pending", ... }
#     }
#   }
#
# All operations are JSON-safe via jq; we degrade gracefully if jq is
# missing (state still tracked, just less rich).
###############################################################################

# ---------------------------------------------------------------------------
# Path helpers.
# ---------------------------------------------------------------------------
ngs_state_path() {
  local sample="${1:?sample required}"
  printf '%s/%s/.state.json' "${work_dir:?work_dir not set}" "$sample"
}

# Has jq? Sample state is best-effort if not.
ngs_has_jq() { command -v jq >/dev/null 2>&1; }

# ---------------------------------------------------------------------------
# Initialise a sample's state file. Idempotent — if the file exists it's
# left alone except for refreshing the layout/group fields.
# ---------------------------------------------------------------------------
ngs_state_init() {
  local sample="${1:?sample required}"
  local file
  file="$(ngs_state_path "$sample")"
  mkdir -p "$(dirname "$file")"

  if [[ ! -f "$file" ]]; then
    if ngs_has_jq; then
      jq -n \
        --arg sid "$sample" \
        --arg ts "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
        --arg layout "${Layout_dict[$sample]:-}" \
        --arg group "${Group_dict[$sample]:-}" \
        --arg batch "${Batch_dict[$sample]:-}" \
        '{sample_id:$sid, created_at:$ts, layout:$layout, group:$group, batch:$batch, stages:{}}' \
        > "$file"
    else
      printf '{"sample_id":"%s","created_at":"%s","layout":"%s","group":"%s","batch":"%s","stages":{}}\n' \
        "$sample" "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
        "${Layout_dict[$sample]:-}" \
        "${Group_dict[$sample]:-}" \
        "${Batch_dict[$sample]:-}" \
        > "$file"
    fi
  fi
}

# ---------------------------------------------------------------------------
# Get a stage's status. Returns one of: missing | pending | running |
# completed | failed.
# ---------------------------------------------------------------------------
ngs_state_stage_status() {
  local sample="${1:?sample required}"
  local stage="${2:?stage required}"
  local file
  file="$(ngs_state_path "$sample")"

  [[ ! -f "$file" ]] && { echo missing; return; }
  if ngs_has_jq; then
    local s
    s="$(jq -r --arg k "$stage" '.stages[$k].status // "missing"' "$file" 2>/dev/null)"
    echo "${s:-missing}"
  else
    if grep -q "\"$stage\"" "$file" 2>/dev/null; then
      echo "completed"   # best-effort fallback; jq is the source of truth
    else
      echo missing
    fi
  fi
}

# Is a stage complete for a sample?
ngs_state_is_complete() {
  [[ "$(ngs_state_stage_status "$1" "$2")" == "completed" ]]
}

# ---------------------------------------------------------------------------
# Begin a stage. Records start time + status="running" so partial failures
# are visible in `ngsmodule status` even before the stage completes.
# ---------------------------------------------------------------------------
ngs_state_stage_begin() {
  local sample="${1:?sample required}"
  local stage="${2:?stage required}"
  local version="${3:-1.0.0}"
  ngs_state_init "$sample"
  ngs_has_jq || return 0
  local file
  file="$(ngs_state_path "$sample")"
  local tmp
  tmp="$(mktemp "$file.XXXXXX")"
  jq --arg k "$stage" --arg ts "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" --arg v "$version" \
     '.stages[$k] = {status:"running", started_at:$ts, version:$v}' \
     "$file" > "$tmp" && mv "$tmp" "$file"
}

# ---------------------------------------------------------------------------
# Mark a stage finished. `outcome` is "completed" or "failed".
#
# Optional rich metadata:
#   --params <json>       # {"min_quality":20, ...}
#   --tools  <json>       # [{"name":"fastp","version":"0.23.4"}, ...]
#   --inputs <json>       # [{"path":"...","sha256":"..."}]
#   --outputs <json>      # [{"path":"...","sha256":"...","kind":"..."}]
#   --error <text>        # only when outcome=failed
# ---------------------------------------------------------------------------
ngs_state_stage_end() {
  local sample="${1:?sample required}"; shift
  local stage="${1:?stage required}"; shift
  local outcome="${1:?outcome required}"; shift

  ngs_has_jq || return 0
  local file
  file="$(ngs_state_path "$sample")"
  [[ ! -f "$file" ]] && return 0

  local params="null" tools="[]" inputs="[]" outputs="[]" err=""
  while (( $# > 0 )); do
    case "$1" in
      --params)  params="$2"; shift 2 ;;
      --tools)   tools="$2"; shift 2 ;;
      --inputs)  inputs="$2"; shift 2 ;;
      --outputs) outputs="$2"; shift 2 ;;
      --error)   err="$2"; shift 2 ;;
      *)         shift ;;
    esac
  done

  local tmp
  tmp="$(mktemp "$file.XXXXXX")"
  jq --arg k "$stage" \
     --arg outcome "$outcome" \
     --arg ts "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
     --argjson params "$params" \
     --argjson tools "$tools" \
     --argjson inputs "$inputs" \
     --argjson outputs "$outputs" \
     --arg err "$err" \
     '
     .stages[$k] = (.stages[$k] // {})
       | .stages[$k].status = $outcome
       | .stages[$k].completed_at = $ts
       | (if .stages[$k].started_at then
            .stages[$k].elapsed_seconds = (
              ($ts | fromdateiso8601) - (.stages[$k].started_at | fromdateiso8601)
            )
          else . end)
       | .stages[$k].params = $params
       | .stages[$k].tools = $tools
       | .stages[$k].inputs = $inputs
       | .stages[$k].outputs = $outputs
       | (if $err != "" then .stages[$k].error = $err else . end)
     ' "$file" > "$tmp" && mv "$tmp" "$file"
}

# ---------------------------------------------------------------------------
# Compute SHA256 for a single file (returns empty string if file missing
# or sha256sum unavailable). Cheap enough that we hash every input/output;
# can be disabled with NGS_NO_HASH=1 for very large cohorts.
# ---------------------------------------------------------------------------
ngs_sha256() {
  local f="$1"
  [[ "${NGS_NO_HASH:-0}" == "1" ]] && { echo ""; return; }
  [[ ! -f "$f" ]] && { echo ""; return; }
  command -v sha256sum >/dev/null 2>&1 || { echo ""; return; }
  sha256sum "$f" 2>/dev/null | awk '{print $1}'
}

# ---------------------------------------------------------------------------
# Build a JSON array of {path, sha256[, kind]} entries. Useful for the
# inputs/outputs columns of stage records.
#
# Usage:
#   ngs_state_files_json /path/a.fq.gz /path/b.fq.gz
#   ngs_state_files_json --kind=trimmed_fastq /path/out_1.fq.gz /path/out_2.fq.gz
# ---------------------------------------------------------------------------
ngs_state_files_json() {
  local kind=""
  if [[ "$1" == --kind=* ]]; then
    kind="${1#--kind=}"
    shift
  fi
  local out="["
  local first=1
  for f in "$@"; do
    [[ -z "$f" ]] && continue
    local sha
    sha="$(ngs_sha256 "$f")"
    if (( ! first )); then out+=","; fi
    out+=$(printf '{"path":"%s","sha256":"%s"' "$f" "$sha")
    [[ -n "$kind" ]] && out+=$(printf ',"kind":"%s"' "$kind")
    out+="}"
    first=0
  done
  out+="]"
  printf '%s' "$out"
}

# ---------------------------------------------------------------------------
# Record a tool's resolved path + version in JSON form. Used as input to
# ngs_state_stage_end --tools.
#
# Usage:
#   tools_json="$(ngs_state_tools_json fastp fastqc samtools)"
# ---------------------------------------------------------------------------
ngs_state_tools_json() {
  local out="["
  local first=1
  for tool in "$@"; do
    local path version
    path="$(command -v "$tool" 2>/dev/null || echo "")"
    if [[ -n "$path" ]]; then
      version="$("$tool" --version 2>&1 | head -1 | awk '{print $NF}')"
    else
      version=""
    fi
    if (( ! first )); then out+=","; fi
    out+=$(printf '{"name":"%s","path":"%s","version":"%s"}' "$tool" "$path" "$version")
    first=0
  done
  out+="]"
  printf '%s' "$out"
}

# ---------------------------------------------------------------------------
# Emit a Provenance Receipt (S3): a self-contained JSON document for the
# stage that lets a future reader reproduce the run.
#
# Receipts are written next to the state file at:
#   $work_dir/$sample/receipts/${stage}_${timestamp}.json
#
# They are append-only; running a stage twice produces two receipts, so
# you can audit how a result was reached even after force-rerun.
# ---------------------------------------------------------------------------
ngs_state_emit_receipt() {
  local sample="${1:?sample required}"
  local stage="${2:?stage required}"
  ngs_has_jq || return 0

  local file
  file="$(ngs_state_path "$sample")"
  [[ ! -f "$file" ]] && return 0

  local sample_dir
  sample_dir="$(dirname "$file")"
  local rdir="$sample_dir/receipts"
  mkdir -p "$rdir"

  local ts
  ts="$(date -u +%Y%m%dT%H%M%SZ)"
  local out="$rdir/${stage}_${ts}.json"

  jq --arg sample "$sample" \
     --arg stage "$stage" \
     --arg generated_at "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" \
     --arg host "$(hostname 2>/dev/null || echo unknown)" \
     --arg user "$(id -un 2>/dev/null || echo unknown)" \
     --arg git "$(git -C "${NGSMODULE_ROOT:-/}" rev-parse --short HEAD 2>/dev/null || echo unknown)" \
     '{
       sample: $sample,
       stage: $stage,
       generated_at: $generated_at,
       host: $host,
       user: $user,
       ngsmodule_git: $git,
       state: .stages[$stage]
     }' "$file" > "$out" 2>/dev/null

  emit_artifact provenance_receipt "$out" 2>/dev/null || true
}
