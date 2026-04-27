# Pipeline Migration Guide

A practical walkthrough for porting a legacy script under `GeneralSteps/`,
`Analysis/`, or `SingleCellPipe/` to the new framework under
`pipelines/core/<Pipeline>/`. Every step here was used to migrate
`postAlignmentQC`, `Quantification`, and `VariantCalling`.

## What you get for free after migrating

- Auto-discovery: `ngsmodule list` shows the pipeline.
- Schema validation: `ngsmodule lint` and the orchestrator both check params.
- Resume / skip / force: handled by `.state.json`, no per-pipeline lockfile.
- Cluster execution: same code runs locally, on Slurm, SGE, or LSF via `--profile`.
- Benchmarks: `/usr/bin/time -v` captured automatically when present.
- HTML report: status grid, gantt, resource peaks, tool tally, failures.
- CI test fixture: `ngsmodule test <pipeline>` runs in `--dry-run` mode.

## Anatomy of a migrated pipeline

```
pipelines/core/<Pipeline>/
├── meta.yml          # display name, version, requires:, schema, resources
├── env.yml           # conda environment (versions pinned)
├── pipeline.sh       # the orchestrator-driven script
└── tests/
    ├── test_config.sh
    ├── test_samples.csv
    └── test_data/.gitkeep
```

## Step-by-step

### 1. Identify the bash script's contract

Find the legacy `*.sh` and answer:

| Question                                          | Where to look                                  |
| ------------------------------------------------- | ---------------------------------------------- |
| What inputs does it consume?                      | `bam=$dir/...` style assignments               |
| What outputs does it produce?                     | `mkdir -p` / `>` redirects / `-o` flags        |
| Which binaries does it invoke?                    | `command_name &>/dev/null` precheck blocks     |
| Which env vars (parameters) does it read?         | `${Aligner}`, `${SequenceType}`, `${threads}`  |
| What pipelines must run first?                    | "Cannot find ... Please run X first" errors    |

Write these on paper before touching any framework files.

### 2. Create the directory + meta.yml

```bash
mkdir -p pipelines/core/MyPipe/tests/test_data
touch    pipelines/core/MyPipe/tests/test_data/.gitkeep
```

`meta.yml` is the framework's source of truth. Use the schema from a
sibling pipeline (`postAlignmentQC/meta.yml` is a good template) and fill in:

```yaml
name: MyPipe
display_name: Friendly name
category: QC | Alignment | Quantification | VariantCalling | Analysis
version: 1.0.0
description: >
  One paragraph. The web UI shows this verbatim.
requires:
  - Alignment            # auto-resolved per sample, omit for root pipelines
tools:
  - some-binary >= 1.0   # informational only — used by `ngsmodule info`
default_params:
  Aligner: STAR
  threads: 4
params_schema:
  Aligner:
    type: enum
    values: [STAR, bwa-mem2, hisat2]
    required: true
    default: STAR
    description: Read aligner used upstream.
  threads:
    type: int
    min: 1
    max: 64
    default: 4
resources:                # consumed by --profile slurm/sge/lsf
  mem: 16G
  time: "2:00:00"
  cpus_per_task: 4
  partition: short
inputs:
  kind: aligned_bam
  pattern: "{sample}.{aligner}.sorted.bam"
outputs:
  - kind: my_output
    pattern: "{sample}.something.txt"
```

Schema types supported by `pipelines/lib/schema.sh`:
`enum`, `int`, `float`, `bool`, `path`, `string` (with optional `pattern`).
For `int` / `float` you can set `min` / `max`.

### 3. Add tool modules (if needed)

Modules live under `pipelines/modules/<tool>.sh`. Each module is a small
bash file exposing `module_<tool>_<action>` functions with **named arguments**:

```bash
#!/usr/bin/env bash
source "$(dirname "${BASH_SOURCE[0]}")/_lib.sh"

module_mytool_run() {
  local in="" out="" threads="${threads:-4}"
  while (( $# > 0 )); do
    case "$1" in
      --in)      in="$2"; shift 2 ;;
      --out)     out="$2"; shift 2 ;;
      --threads) threads="$2"; shift 2 ;;
      *) echo "module_mytool_run: unknown arg $1" >&2; return 2 ;;
    esac
  done
  : "${in:?--in required}"; : "${out:?--out required}"

  module_run "mytool" \
    mytool --threads "$threads" -i "$in" -o "$out"
}
```

`module_run` does three things:
1. Tracks the tool name in `_MODULES_TOOLS_USED` for the state.json report.
2. Honours `--dry-run` (prints `::dry-run::` instead of executing).
3. Wraps the command with `/usr/bin/time -v` when benchmarks are enabled.

If the underlying CLI is multi-command (e.g. `bcftools mpileup | call`),
use `bash -c "..."` inside `module_run` and set
`NGS_MODULE_CURRENT_TOOL=bcftools` so the tool name is captured correctly.

### 4. Write `pipeline.sh`

The contract is a single function `run_<Pipeline>_for_sample <sample>`:

```bash
#!/usr/bin/env bash
set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_NGS_LIB_DIR="$(cd "$_PIPELINE_DIR/../../lib" && pwd)"
export _NGS_LIB_DIR

source "$_NGS_LIB_DIR/orchestrator.sh"

_NGS_MODULES_DIR="$(cd "$_PIPELINE_DIR/../../modules" && pwd)"
source "$_NGS_MODULES_DIR/_lib.sh"
module_load mytool   # one per binary you'll invoke

run_MyPipe_for_sample() {
  local sample="${1:?sample required}"
  local sample_dir="$work_dir/$sample"
  module_reset                                  # fresh per-sample tool tally

  # -- Locate inputs --------------------------------------------------
  local in_file="$sample_dir/Alignment-${Aligner}/${sample}.${Aligner}.sorted.bam"
  if [[ ! -f "$in_file" ]] && ! is_dry_run; then
    log_error "[$sample] missing input: $in_file"
    return 1
  fi

  # -- Run the work via modules --------------------------------------
  emit_status "[$sample] MyPipe: doing the thing"
  local out_file="$sample_dir/MyPipe/${sample}.out"
  mkdir -p "$(dirname "$out_file")"
  module_mytool_run --in "$in_file" --out "$out_file" --threads "${threads:-4}"

  # -- Record provenance + artifacts ---------------------------------
  if ! is_dry_run; then
    emit_artifact my_output "$out_file"
  fi

  local tools_json inputs_json outputs_json params_json
  tools_json="$(module_tools_json)"
  inputs_json="$(ngs_state_files_json --kind=aligned_bam "$in_file")"
  outputs_json="$(ngs_state_files_json --kind=my_output "$out_file")"
  params_json="$(printf '{"Aligner":"%s","threads":%s}' "${Aligner:-STAR}" "${threads:-4}")"

  ngs_state_stage_end "$sample" MyPipe completed \
    --params "$params_json" --tools "$tools_json" \
    --inputs "$inputs_json" --outputs "$outputs_json"
}

# Standalone-execution shim — keep this block; it lets users run
# `pipelines/core/MyPipe/pipeline.sh` directly with -c ConfigFile.
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  if [[ -n "${ConfigFile:-}" ]] && [[ -f "${ConfigFile}" ]]; then
    source "$ConfigFile"
  fi
  if [[ -n "${SampleInfoFile:-}" ]] && [[ -f "${SampleInfoFile}" ]]; then
    ngs_load_sample_info "$SampleInfoFile"
  fi
  if [[ -n "${rawdata_dir:-}" ]] && [[ -n "${work_dir:-}" ]]; then
    if [[ ! -d "$work_dir" ]] || [[ -z "$(ls -A "$work_dir" 2>/dev/null)" ]]; then
      ngs_build_workdir "$rawdata_dir" "$work_dir"
    fi
  fi
  pipeline_run MyPipe
fi
```

What the orchestrator does for you (you don't write any of this):
- Iterates over `samples[]` with `ntask_per_run` concurrency
- Schema-validates env vars before the loop
- Auto-fills defaults from `params_schema:`
- Resolves `requires:` per sample (skipping prereqs that already completed)
- Emits `::status::` / `::progress::` / `::artifact::` events for the web UI
- Honors `--retry`, `--force`, `--no-deps`, `--dry-run`
- Captures wall/CPU/RSS via `/usr/bin/time -v` when present
- Persists everything to `$work_dir/<sample>/.state.json`

### 5. Add the test fixture

```bash
# pipelines/core/MyPipe/tests/test_samples.csv
RunID,SampleID,Layout,Group
runA,test_sample_PE,PE,test
```

```bash
# pipelines/core/MyPipe/tests/test_config.sh
rawdata_dir="${TEST_DIR:-.}/test_data"
work_dir="${TEST_WORK:-/tmp/ngsmodule-test}/MyPipe"
SampleInfoFile="${TEST_DIR:-.}/test_samples.csv"

ntask_per_run=1
threads=2
Aligner=STAR
```

Empty `tests/test_data/.gitkeep` and you're done — `ngsmodule test MyPipe`
will run the pipeline in dry-run mode and pass if `pipeline.sh` doesn't error.

### 6. Validate

```bash
./ngsmodule lint MyPipe          # offline schema + meta.yml check
./ngsmodule test MyPipe          # dry-run regression
./ngsmodule info MyPipe          # render meta.yml
./ngsmodule list                 # confirm the pipeline shows up
```

## Common gotchas

- **`requires:` is per-sample, not global.** If a sample's prereq is already
  in `.state.json` as `completed`, it gets skipped automatically.
- **Don't `cd` inside `run_*_for_sample`.** Return to the previous dir before
  the function exits, or — preferred — pass absolute paths to module functions.
- **Quote everything that flows through heredocs.** Profile scripts (`slurm.sh`,
  `sge.sh`, `lsf.sh`) generate batch scripts via `cat <<EOF`; values with spaces
  or single quotes will break unquoted expansions.
- **Tool name with hyphens.** `module_load bwa-mem2` won't work; the file is
  `pipelines/modules/bwa_mem2.sh` and the load guard converts hyphens to
  underscores. Source-side code uses `_MODULES_TOOLS_USED[bwa-mem2]=1` (key
  remains hyphenated) for the report tally.
- **Schema enum quoting.** YAML inline lists strip surrounding quotes per
  element: `values: [STAR, "bwa-mem2", hisat2]` is parsed as
  `STAR,bwa-mem2,hisat2`. Pin known-tricky values explicitly.

## Reference: existing migrations

| Pipeline         | Modules used                                           |
| ---------------- | ------------------------------------------------------ |
| preAlignmentQC   | `fastp`, `fastqc`                                      |
| Alignment        | `star`, `bwa_mem2`, `hisat2`, `samtools`               |
| postAlignmentQC  | `rseqc`, `preseq`, `mosdepth`, `goleft`                |
| Quantification   | `featurecounts`                                        |
| VariantCalling   | `bcftools`, `tabix`                                    |

Each is a self-contained example of the pattern above. Read the `pipeline.sh`
of whichever is closest to what you're migrating before writing your own.
