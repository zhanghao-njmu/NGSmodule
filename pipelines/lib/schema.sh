#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/schema.sh
#
# Pipeline parameter schema validation, inspired by nf-core's
# nextflow_schema.json but kept YAML-native and pure-bash.
#
# Schema lives inline in meta.yml under `params_schema:` and is consumed
# by:
#   - pipeline_validate_params      (runtime — abort early on bad config)
#   - ngsmodule lint                (offline — CI/dev quality gate)
#   - backend/init_pipeline_templates.py (auto-generates web UI form widgets)
#
# Schema YAML grammar (intentionally flat — one parser, no PyYAML):
#
#   params_schema:
#     Aligner:
#       type: enum
#       values: [STAR, bwa-mem2, hisat2, bismark]
#       default: STAR
#       required: true
#       description: Read aligner to use.
#     threads:
#       type: int
#       min: 1
#       max: 256
#       default: 4
#     min_quality:
#       type: int
#       min: 0
#       max: 50
#     SequenceType:
#       type: enum
#       values: [rna, dna]
#     output_prefix:
#       type: string
#       pattern: '^[A-Za-z0-9_]+$'
#
# Supported types: enum, int, float, bool, string, path
# Supported attrs: type, values, min, max, default, required, pattern, description
###############################################################################

# ---------------------------------------------------------------------------
# Parse params_schema from meta.yml into an associative array.
#
# Output keys are "<param>.<attr>" so the caller can introspect.
#
# Usage:
#   declare -A SCHEMA
#   pipeline_schema_load path/to/meta.yml SCHEMA
#   echo "${SCHEMA[Aligner.type]}"      # → enum
#   echo "${SCHEMA[Aligner.values]}"    # → STAR,bwa-mem2,hisat2,bismark
#   echo "${SCHEMA[__params__]}"        # → "Aligner threads min_quality ..."
# ---------------------------------------------------------------------------
pipeline_schema_load() {
  local file="${1:?meta yml path required}"
  local -n _out="$2"
  [[ ! -f "$file" ]] && return 0

  local in_block=0 in_param=0 param="" line key val params_list=""
  while IFS= read -r line || [[ -n "$line" ]]; do
    # Strip trailing CR (Windows line endings).
    line="${line%$'\r'}"
    # Skip pure comment lines (but preserve indentation for blank-line context).
    [[ "$line" =~ ^[[:space:]]*# ]] && continue

    # Block boundary detection.
    if [[ "$line" =~ ^params_schema:[[:space:]]*$ ]]; then
      in_block=1; in_param=0; param=""; continue
    fi
    if (( in_block )) && [[ "$line" =~ ^[A-Za-z_] ]]; then
      # Left-aligned new top-level key → leaving the schema block.
      in_block=0; in_param=0; param=""; continue
    fi
    (( in_block )) || continue

    # Param header: 2-space indent, "<name>:" with nothing after.
    if [[ "$line" =~ ^[[:space:]]{2}([A-Za-z_][A-Za-z_0-9]*):[[:space:]]*$ ]]; then
      param="${BASH_REMATCH[1]}"
      in_param=1
      params_list+="$param "
      continue
    fi

    # Param attribute: 4+ space indent, "<key>: <value>".
    if (( in_param )) && [[ "$line" =~ ^[[:space:]]{4,}([A-Za-z_][A-Za-z_0-9]*):[[:space:]]*(.*)$ ]]; then
      key="${BASH_REMATCH[1]}"
      val="${BASH_REMATCH[2]}"
      # Strip surrounding quotes.
      val="${val%\"}"; val="${val#\"}"
      val="${val%\'}"; val="${val#\'}"
      # Strip surrounding [ ] for inline lists, normalise comma+space → comma,
      # and strip per-element quotes so `[STAR, "bwa-mem2", hisat2]` becomes
      # `STAR,bwa-mem2,hisat2`.
      if [[ "$val" =~ ^\[(.*)\]$ ]]; then
        val="${BASH_REMATCH[1]}"
        val="${val// /}"
        val="${val//\"/}"
        val="${val//\'/}"
      fi
      _out["$param.$key"]="$val"
      continue
    fi
  done < "$file"

  _out["__params__"]="${params_list% }"
}

# ---------------------------------------------------------------------------
# Validate a single param value against its schema entry. Echoes a
# human-readable error to stderr and returns non-zero on failure.
# ---------------------------------------------------------------------------
_schema_check_one() {
  local -n _s="$1"
  local param="$2"
  local value="$3"

  local type="${_s[$param.type]:-string}"
  case "$type" in
    enum)
      local values="${_s[$param.values]:-}"
      [[ -z "$values" ]] && return 0
      local IFS=','
      local v ok=0
      for v in $values; do
        [[ "$value" == "$v" ]] && { ok=1; break; }
      done
      if (( ! ok )); then
        printf '  %s=%s — not in allowed values [%s]\n' "$param" "$value" "$values" >&2
        return 1
      fi
      ;;
    int)
      if ! [[ "$value" =~ ^-?[0-9]+$ ]]; then
        printf '  %s=%s — not an integer\n' "$param" "$value" >&2
        return 1
      fi
      local min="${_s[$param.min]:-}" max="${_s[$param.max]:-}"
      if [[ -n "$min" ]] && (( value < min )); then
        printf '  %s=%s — below min (%s)\n' "$param" "$value" "$min" >&2
        return 1
      fi
      if [[ -n "$max" ]] && (( value > max )); then
        printf '  %s=%s — above max (%s)\n' "$param" "$value" "$max" >&2
        return 1
      fi
      ;;
    float)
      if ! [[ "$value" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
        printf '  %s=%s — not a number\n' "$param" "$value" >&2
        return 1
      fi
      ;;
    bool)
      case "${value,,}" in
        true|false|yes|no|0|1) ;;
        *)
          printf '  %s=%s — not a boolean\n' "$param" "$value" >&2
          return 1
          ;;
      esac
      ;;
    path)
      # Path may not exist yet (it could be an output) — only validate
      # syntactic shape: not empty, not a tab-or-newline (those would
      # break shell-escaping downstream). Bash strings can't actually
      # contain NULs (the previous null check matched every value
      # because bash strips $'\0' from the glob pattern, leaving `*`).
      if [[ -z "$value" ]] || [[ "$value" == *$'\t'* ]] || [[ "$value" == *$'\n'* ]]; then
        printf '  %s=%s — not a valid path\n' "$param" "$value" >&2
        return 1
      fi
      ;;
    string)
      local pattern="${_s[$param.pattern]:-}"
      if [[ -n "$pattern" ]] && ! [[ "$value" =~ $pattern ]]; then
        printf '  %s=%s — does not match pattern /%s/\n' "$param" "$value" "$pattern" >&2
        return 1
      fi
      ;;
    *)
      printf '  %s — unknown schema type "%s" (skipping)\n' "$param" "$type" >&2
      ;;
  esac
  return 0
}

# ---------------------------------------------------------------------------
# Validate the running environment's variables against a pipeline's
# params_schema. Aborts (returns 1) on the first batch of failures so
# users see *all* their config errors at once, not one-at-a-time.
#
# Usage:
#   pipeline_validate_params Alignment   # uses $PIPELINES_DIR/Alignment/meta.yml
# ---------------------------------------------------------------------------
pipeline_validate_params() {
  local pipeline="${1:?pipeline name required}"
  local meta="${PIPELINES_DIR:-$NGSMODULE_ROOT/pipelines/core}/$pipeline/meta.yml"
  [[ ! -f "$meta" ]] && return 0

  declare -A SCHEMA
  pipeline_schema_load "$meta" SCHEMA

  local params="${SCHEMA[__params__]:-}"
  [[ -z "$params" ]] && return 0

  local errors=0 missing=0 p value required default
  for p in $params; do
    required="${SCHEMA[$p.required]:-false}"
    default="${SCHEMA[$p.default]:-}"
    value="${!p:-}"

    # If unset, fall back to default before deciding "missing".
    if [[ -z "$value" ]] && [[ -n "$default" ]]; then
      printf -v "$p" '%s' "$default"
      export "${p?}"
      value="$default"
    fi

    if [[ -z "$value" ]]; then
      case "${required,,}" in
        true|yes|1)
          printf '  %s — required but unset\n' "$p" >&2
          missing=$((missing + 1))
          ;;
      esac
      continue
    fi

    if ! _schema_check_one SCHEMA "$p" "$value"; then
      errors=$((errors + 1))
    fi
  done

  if (( errors > 0 || missing > 0 )); then
    printf '\n[%s] schema validation failed: %d missing, %d invalid\n' \
      "$pipeline" "$missing" "$errors" >&2
    printf 'Edit your config file or set the variables before re-running.\n' >&2
    return 1
  fi
  return 0
}

# ---------------------------------------------------------------------------
# Lint a pipeline's meta.yml for completeness and obvious problems.
# Returns the number of issues found (so you can chain in CI).
#
# Checks:
#   - Required top-level keys: name, version, description
#   - Has env.yml or container declared
#   - default_params keys are reflected in params_schema (warn only)
#   - Schema entries have a valid type
#   - enum schemas declare values
#   - int/float schemas with both min and max satisfy min < max
# ---------------------------------------------------------------------------
pipeline_lint() {
  local pipeline="${1:?pipeline required}"
  local pdir="${PIPELINES_DIR:-$NGSMODULE_ROOT/pipelines/core}/$pipeline"
  local meta="$pdir/meta.yml"
  local issues=0

  if [[ ! -f "$meta" ]]; then
    printf '✗ %s: meta.yml missing\n' "$pipeline" >&2
    return 1
  fi

  # Required scalar keys.
  local k
  for k in name version description; do
    if ! grep -qE "^${k}:" "$meta"; then
      printf '✗ %s: meta.yml missing required key: %s\n' "$pipeline" "$k" >&2
      issues=$((issues + 1))
    fi
  done

  # env.yml or container.
  if [[ ! -f "$pdir/env.yml" ]] && ! grep -qE '^container:' "$meta"; then
    printf '⚠ %s: neither env.yml nor `container:` declared in meta.yml\n' "$pipeline" >&2
    issues=$((issues + 1))
  fi

  # Schema sanity.
  declare -A SCHEMA
  pipeline_schema_load "$meta" SCHEMA
  local params="${SCHEMA[__params__]:-}"
  # Warn if the block is declared but empty (likely an editing accident).
  if grep -qE '^params_schema:[[:space:]]*$' "$meta" && [[ -z "$params" ]]; then
    printf '⚠ %s: params_schema: block declared but contains no parameters\n' \
      "$pipeline" >&2
    issues=$((issues + 1))
  fi
  local p type values minv maxv
  for p in $params; do
    type="${SCHEMA[$p.type]:-string}"
    case "$type" in
      enum|int|float|bool|string|path) ;;
      *)
        printf '✗ %s: param %s has unknown type "%s"\n' "$pipeline" "$p" "$type" >&2
        issues=$((issues + 1))
        ;;
    esac
    if [[ "$type" == "enum" ]]; then
      values="${SCHEMA[$p.values]:-}"
      if [[ -z "$values" ]]; then
        printf '✗ %s: enum param %s missing values:\n' "$pipeline" "$p" >&2
        issues=$((issues + 1))
      fi
    fi
    if [[ "$type" == "int" || "$type" == "float" ]]; then
      minv="${SCHEMA[$p.min]:-}"; maxv="${SCHEMA[$p.max]:-}"
      if [[ -n "$minv" && -n "$maxv" ]]; then
        # Float-safe comparison via awk.
        if ! awk -v a="$minv" -v b="$maxv" 'BEGIN{exit !(a < b)}'; then
          printf '✗ %s: param %s has min (%s) >= max (%s)\n' \
            "$pipeline" "$p" "$minv" "$maxv" >&2
          issues=$((issues + 1))
        fi
      fi
    fi
  done

  # Pipeline.sh exists.
  if [[ ! -f "$pdir/pipeline.sh" ]]; then
    printf '✗ %s: pipeline.sh missing\n' "$pipeline" >&2
    issues=$((issues + 1))
  fi

  # `requires:` entries must point to real pipelines.
  local req prereq_dir
  for req in $(awk '
    /^requires:/        {in_block=1; next}
    /^[A-Za-z_]+:/      {in_block=0}
    in_block && /^[[:space:]]*-/ {
      sub(/^[[:space:]]*-[[:space:]]*/, "")
      print
    }
  ' "$meta"); do
    prereq_dir="${PIPELINES_DIR:-$NGSMODULE_ROOT/pipelines/core}/$req"
    if [[ ! -d "$prereq_dir" ]]; then
      printf '✗ %s: requires unknown pipeline "%s"\n' "$pipeline" "$req" >&2
      issues=$((issues + 1))
    fi
  done

  if (( issues == 0 )); then
    printf '✓ %s: clean\n' "$pipeline"
  else
    printf '✗ %s: %d issue(s)\n' "$pipeline" "$issues" >&2
  fi
  return $issues
}
