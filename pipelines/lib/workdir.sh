#!/usr/bin/env bash
# shellcheck shell=bash
###############################################################################
# pipelines/lib/workdir.sh
#
# Raw-data → work-directory layout management via softlinks.
#
# Design tenet: original FASTQ files are NEVER touched. The work directory
# is a derived space holding symlinks (and later, processed BAM/VCF outputs)
# that can be wiped and rebuilt without touching primary data.
#
# Layout produced:
#
#   $rawdata_dir/                            ← read-only primary data
#     GSE12345_R1_L001_R1_001.fastq.gz
#     GSE12345_R1_L001_R2_001.fastq.gz
#     GSE12345_R2_L001_R1_001.fastq.gz       ← second sequencing run, same sample
#     GSE12345_R2_L001_R2_001.fastq.gz
#
#   $work_dir/                                ← derived; safe to rm -rf
#     Tumor_S01/                              ← SampleID from SampleInfoFile
#       run1_Tumor_S01_1.fq.gz  → ../../$rawdata_dir/GSE12345_R1_L001_R1_001.fastq.gz
#       run1_Tumor_S01_2.fq.gz  → ../../$rawdata_dir/GSE12345_R1_L001_R2_001.fastq.gz
#       run2_Tumor_S01_1.fq.gz  → ../../$rawdata_dir/GSE12345_R2_L001_R1_001.fastq.gz
#       run2_Tumor_S01_2.fq.gz  → ../../$rawdata_dir/GSE12345_R2_L001_R2_001.fastq.gz
#       PreAlignmentQC/, Alignment-bwa/, ...  ← created by pipelines
#
# Key conventions:
#   - run{N}_ prefix supports multi-run samples (auto-merged before mapping)
#   - _1 / _2 suffix for paired-end; bare for single-end
#   - SampleID never derived from filenames; always from SampleInfoFile mapping
#
# Public API:
#   ngs_load_sample_info <csv_path>          → fills Sample_dict, Layout_dict,
#                                              Group_dict, Batch_dict, samples[]
#   ngs_build_workdir <rawdata> <workdir>    → creates symlinks per convention
#   ngs_archive_workdir <workdir>            → atomic rename to bk_<ts>_work
#   ngs_sample_inputs <sample_workdir>       → echoes input fq files for a sample
###############################################################################

# ---------------------------------------------------------------------------
# Suffix patterns. Override before calling ngs_build_workdir if your raw
# files use a different convention. Defaults match Illumina BaseSpace.
# ---------------------------------------------------------------------------
: "${R1_SufixPattern:=_1.fq.gz|_R1.fq.gz|_R1_001.fastq.gz|_1.fastq.gz}"
: "${R2_SufixPattern:=_2.fq.gz|_R2.fq.gz|_R2_001.fastq.gz|_2.fastq.gz}"
: "${SE_SufixPattern:=.fq.gz|.fastq.gz}"
: "${RunIDPattern:=[A-Za-z0-9_-]+}"

# Globals populated by ngs_load_sample_info (associative arrays).
declare -gA Sample_dict      # RunID → SampleID
declare -gA Layout_dict      # SampleID → PE | SE
declare -gA Group_dict       # SampleID → group label (e.g. Treatment / Control)
declare -gA Batch_dict       # SampleID → batch label
declare -ga samples          # ordered unique sample list

# ---------------------------------------------------------------------------
# Parse the SampleInfoFile (CSV with header).
#
# Required columns: RunID, SampleID, Layout
# Optional columns: Group, Batch (or BatchID)
#
# Lines starting with # are skipped. The header row is auto-detected by
# matching "RunID" case-insensitively in the first column.
# ---------------------------------------------------------------------------
ngs_load_sample_info() {
  local csv="${1:?ngs_load_sample_info: CSV path required}"
  if [[ ! -f "$csv" ]]; then
    echo "ERROR: SampleInfoFile not found: $csv" >&2
    return 1
  fi

  Sample_dict=()
  Layout_dict=()
  Group_dict=()
  Batch_dict=()
  samples=()

  local header_seen=0
  local idx_run idx_sample idx_layout idx_group idx_batch
  local line
  while IFS= read -r line || [[ -n "$line" ]]; do
    # strip CR
    line="${line%$'\r'}"
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^# ]] && continue

    if (( ! header_seen )); then
      # detect header
      local first="${line%%,*}"
      if [[ "${first,,}" == *"runid"* ]]; then
        local fields=()
        IFS=',' read -ra fields <<< "$line"
        local i
        for i in "${!fields[@]}"; do
          local h="${fields[$i],,}"
          h="${h// /}"
          case "$h" in
            runid)             idx_run=$i ;;
            sampleid|sample)   idx_sample=$i ;;
            layout)            idx_layout=$i ;;
            group|groupname)   idx_group=$i ;;
            batch|batchid)     idx_batch=$i ;;
          esac
        done
        header_seen=1
        continue
      fi
      # No header — assume canonical order
      idx_run=0; idx_sample=1; idx_layout=2; idx_group=3; idx_batch=4
      header_seen=1
    fi

    local row=()
    IFS=',' read -ra row <<< "$line"
    local run="${row[$idx_run]:-}"
    local samp="${row[$idx_sample]:-}"
    local layout="${row[${idx_layout:-2}]:-}"
    local group="${row[${idx_group:-3}]:-}"
    local batch="${row[${idx_batch:-4}]:-}"

    # trim
    run="${run// /}"; samp="${samp// /}"; layout="${layout// /}"
    group="${group// /}"; batch="${batch// /}"

    if [[ -z "$run" ]] || [[ -z "$samp" ]]; then
      continue
    fi
    if [[ "$layout" != "PE" ]] && [[ "$layout" != "SE" ]]; then
      echo "ERROR: SampleInfoFile row for $run has Layout='$layout' (must be PE or SE)" >&2
      return 1
    fi

    Sample_dict["$run"]="$samp"
    Layout_dict["$samp"]="$layout"
    [[ -n "$group" ]] && Group_dict["$samp"]="$group"
    [[ -n "$batch" ]] && Batch_dict["$samp"]="$batch"

    # append unique to samples[]
    local seen=0 s
    for s in "${samples[@]}"; do [[ "$s" == "$samp" ]] && { seen=1; break; }; done
    (( seen == 0 )) && samples+=("$samp")
  done < "$csv"

  if (( ${#Sample_dict[@]} == 0 )); then
    echo "ERROR: SampleInfoFile produced no usable rows: $csv" >&2
    return 1
  fi
  return 0
}

# ---------------------------------------------------------------------------
# Archive an existing work directory before recreating it.
# Preserves data via rename (atomic) so a failed build doesn't lose
# previous results.
# ---------------------------------------------------------------------------
ngs_archive_workdir() {
  local workdir="${1:?ngs_archive_workdir: workdir required}"
  [[ ! -d "$workdir" ]] && return 0
  local ts
  ts="$(date +%Y%m%d%H%M%S)"
  local parent
  parent="$(dirname "$workdir")"
  local base
  base="$(basename "$workdir")"
  mv "$workdir" "$parent/bk_${ts}_${base}"
  mkdir -p "$workdir"
}

# ---------------------------------------------------------------------------
# Build the work directory by symlinking matched raw files.
#
# Prerequisites: ngs_load_sample_info must have been called.
#
# Args:
#   $1  rawdata_dir (absolute path; must contain FASTQ files)
#   $2  work_dir    (will be created; existing contents are preserved)
# ---------------------------------------------------------------------------
ngs_build_workdir() {
  local rawdata="${1:?rawdata_dir required}"
  local workdir="${2:?work_dir required}"

  if [[ ! -d "$rawdata" ]]; then
    echo "ERROR: rawdata_dir not found: $rawdata" >&2
    return 1
  fi
  if (( ${#Sample_dict[@]} == 0 )); then
    echo "ERROR: ngs_load_sample_info must run first" >&2
    return 1
  fi

  mkdir -p "$workdir"

  local grep_pattern="(${RunIDPattern}(${R1_SufixPattern})\$)|(${RunIDPattern}(${R2_SufixPattern})\$)|(${RunIDPattern}(${SE_SufixPattern})\$)"

  local -a found_files=()
  while IFS= read -r f; do
    found_files+=("$f")
  done < <(find "$rawdata" -type f 2>/dev/null | grep -P "$grep_pattern" | sort)

  if (( ${#found_files[@]} == 0 )); then
    # In dry-run mode the test fixtures intentionally have no real FASTQs;
    # downgrade the hard error to a warning so pipelines can still exercise
    # their schema/module/state flows.
    if declare -f is_dry_run >/dev/null && is_dry_run; then
      echo "WARN: No raw files under $rawdata (dry-run; skipping workdir build)" >&2
      return 0
    fi
    echo "ERROR: No raw files matched the suffix pattern under $rawdata" >&2
    echo "       R1: $R1_SufixPattern" >&2
    echo "       R2: $R2_SufixPattern" >&2
    echo "       SE: $SE_SufixPattern" >&2
    return 1
  fi

  local linked=0 skipped=0
  local file
  for file in "${found_files[@]}"; do
    local fname="${file##*/}"
    local r1_match r2_match se_match
    r1_match=$(echo "$fname" | grep -oP "($R1_SufixPattern)\$" || true)
    r2_match=$(echo "$fname" | grep -oP "($R2_SufixPattern)\$" || true)
    se_match=$(echo "$fname" | grep -oP "($SE_SufixPattern)\$" || true)

    local matched_suffix run_id sample_id fq_layout fq_target
    matched_suffix=""
    for sfx in "$r1_match" "$r2_match" "$se_match"; do
      [[ -z "$sfx" ]] && continue
      local stripped="${fname%$sfx}"
      run_id=$(echo "$stripped" | grep -oP "^${RunIDPattern}\$" || true)
      [[ -z "$run_id" ]] && continue
      if [[ -n "${Sample_dict[$run_id]:-}" ]]; then
        matched_suffix="$sfx"
        break
      fi
    done

    if [[ -z "$matched_suffix" ]] || [[ -z "$run_id" ]]; then
      skipped=$((skipped + 1))
      continue
    fi

    sample_id="${Sample_dict[$run_id]}"
    if [[ "$matched_suffix" == "$r1_match" ]]; then
      fq_target="run1_${sample_id}_1.fq.gz"; fq_layout="PE"
    elif [[ "$matched_suffix" == "$r2_match" ]]; then
      fq_target="run1_${sample_id}_2.fq.gz"; fq_layout="PE"
    else
      fq_target="run1_${sample_id}.fq.gz"; fq_layout="SE"
    fi

    local declared_layout="${Layout_dict[$sample_id]:-}"
    if [[ -n "$declared_layout" ]] && [[ "$declared_layout" != "$fq_layout" ]]; then
      echo "ERROR: $fname is $fq_layout but SampleInfoFile says $sample_id is $declared_layout" >&2
      return 1
    fi

    mkdir -p "$workdir/$sample_id"
    local target="$workdir/$sample_id/$fq_target"
    if [[ ! -e "$target" ]]; then
      ln -s "$file" "$target"
    else
      # multi-run: count existing run* with the same suffix to derive next index
      local suffix_part="${fq_target#run1_}"
      local n
      n=$(find "$workdir/$sample_id" -maxdepth 1 -name "run*_${suffix_part}" 2>/dev/null | wc -l)
      ln -s "$file" "$workdir/$sample_id/run$((n + 1))_${suffix_part}"
    fi
    linked=$((linked + 1))
  done

  echo "ngs_build_workdir: linked=$linked skipped=$skipped samples=${#samples[@]}"
}

# ---------------------------------------------------------------------------
# List the input fastq paths for a sample, in deterministic run order.
# Used by step scripts that need to concatenate multi-run files.
#
# Args:
#   $1  sample workdir (e.g. $work_dir/Tumor_S01)
#   $2  layout (PE | SE)
#   $3  output variable name for read1 array (e.g. r1_files)
#   $4  output variable name for read2 array (PE only)
# ---------------------------------------------------------------------------
ngs_sample_inputs() {
  local sample_dir="${1:?sample_dir required}"
  local layout="${2:?layout required}"
  local -n _r1="$3"
  _r1=()
  if [[ "$layout" == "PE" ]]; then
    local -n _r2="$4"
    _r2=()
    while IFS= read -r f; do _r1+=("$f"); done < <(find "$sample_dir" -maxdepth 1 -name "run*_*_1.fq.gz" | sort)
    while IFS= read -r f; do _r2+=("$f"); done < <(find "$sample_dir" -maxdepth 1 -name "run*_*_2.fq.gz" | sort)
  else
    while IFS= read -r f; do _r1+=("$f"); done < <(find "$sample_dir" -maxdepth 1 -name "run*_*.fq.gz" \
      ! -name "*_1.fq.gz" ! -name "*_2.fq.gz" | sort)
  fi
}
