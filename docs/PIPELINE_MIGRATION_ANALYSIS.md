# Migrating NGSmodule to Snakemake or Nextflow — Tradeoff Analysis

This document evaluates whether to migrate NGSmodule's existing
shell-script pipeline subsystem to a workflow manager (Snakemake or
Nextflow), and what would be lost or gained.

**Verdict (TL;DR): Don't migrate yet.** The current shell scripts
encode lab-specific tribal knowledge that is cheap to fork and modify,
and they integrate trivially with the FastAPI + Celery web layer via
`subprocess.Popen`. A migration would cost 2–4 person-weeks and
introduce a second orchestrator on top of Celery for limited near-term
gain. Prefer **stage A + B optimisations** (this branch) and revisit
migration only once a concrete scaling driver appears (e.g., HPC/cloud
bursting, multi-tenant resource isolation, or a regulatory provenance
requirement).

---

## What we have today

A 514-line `NGSmodule` bash dispatcher backed by ~30 shell scripts under
`GeneralSteps/`, `Analysis/`, `PreparationSteps/`, `PrefetchData/`,
`SingleCellPipe/`. Key non-obvious features:

1. **FIFO-based job semaphore** for sample-level parallelism (no GNU
   parallel, no scheduler).
2. **Logfile-based resume**: completion detected by grepping
   `"NGSmodule finished the job"`; errors detected by regex
   `error|fatal|corrupt|interrupt|EOFException`.
3. **Per-sample retry loop** with configurable `retry` count.
4. **Two-tier configuration**: human-readable `ConfigFile.sh` with ~120
   parameters + `SampleInfoFile.csv` mapping `RunID → SampleID →
   Group → Layout(PE/SE) → BatchID`.
5. **iGenomes layout convention** auto-resolves index paths:
   `$iGenomes_Dir/$Species/$Source/$Build/Sequence/{BWA,STAR,Bismark}Index/`.
6. **Automatic dedup logic** keyed by `SequenceType`
   (`rna`/`dna`/`BSdna`).
7. **Multi-run concatenation**: `run1_, run2_` files for the same
   sample auto-merged before processing.
8. **`trap_add` composable signal handler** with process-group kill on
   Ctrl+C.
9. **Conda env per workflow mode** activated by `CheckENV.sh`.
10. **Color-coded progress bar** in interactive runs.

These are **value-additive features** built up over multiple lab
projects, not generic infrastructure.

---

## What you would gain from Nextflow

| Capability | Nextflow built-in | Current state |
|------------|-------------------|---------------|
| DAG visualization | `-with-dag flow.html` | None |
| HPC / cloud executor | `slurm`, `awsbatch`, `googlebatch`, `k8s`, `local` profiles | Local Linux only |
| Per-process resource decl. | `cpus 16; memory '64.GB'; time '12.h'` | Implicit `threads` cap at 32 |
| Container per process | `container 'biocontainers/fastp:0.23.4'` | Conda env per mode |
| Resume by content hash | `-resume` (Merkle tree) | Logfile string match |
| MultiQC reports | `nf-core/multiqc` module | Manual postcheck |
| Standard wrapper library | `nf-core/modules` (1000+ tools) | Hand-written wrappers |
| Test framework | `nf-test` | None |
| Provenance graph | Auto-generated | None |

## What you would gain from Snakemake

| Capability | Snakemake built-in | Current state |
|------------|--------------------|---------------|
| Python-native rules | `Snakefile` is Python DSL | Bash |
| Automatic input/output dependency | Built-in (rule body) | Explicit chain |
| `--rerun-incomplete` | Built-in | Manual |
| Cluster profiles | `--profile slurm` | Local only |
| Container support | `singularity:` per rule | Conda only |
| Wrapper repo | `snakemake-wrappers` | Hand-written |
| Local Python embedding | `import snakemake` | N/A |

---

## What you would lose

### 1. The "fork-and-tweak" workflow

A wet-lab bioinformatician can today open `preAlignmentQC.sh`, see 445
lines of readable bash, and change `--qualified_quality_phred 20` to
`30` in 30 seconds. With Nextflow they would need to:

- Learn Groovy DSL2 syntax
- Understand `process` / `workflow` / `params` scopes
- Re-run `nextflow run main.nf` with the right `-resume` semantics

For labs with rotating students this **raises the barrier to
contribution** measurably. Not a deal-breaker, but real.

### 2. Implicit defaults

Today `LoadConfig.sh` derives ~30 values from a few inputs (e.g.,
`Deduplication="automatic"` resolves to `True/False` based on
`SequenceType`). Nextflow params are explicit; you would either
re-implement the derivation in a custom Groovy `init` block or push
that complexity onto users who must set every value.

### 3. Index path convention

iGenomes path resolution is encoded in `LoadConfig.sh`:

```bash
genome="$iGenomes_Dir/$Species/$Source/$Build/Sequence/WholeGenomeFasta/genome.fa"
bwa_index="$iGenomes_Dir/$Species/$Source/$Build/Sequence/BWAIndex/genome.fa"
```

Replicating this in Nextflow requires either many `params` or a custom
helper module. Doable but tedious to keep in sync if Illumina ever
restructures the iGenomes tree.

### 4. The web integration becomes more indirect

Today:

```
FastAPI → Celery worker → subprocess.Popen("bash GeneralSteps/X.sh")
                          ↳ tool runs, writes log
                          ↳ worker tails log, parses progress
```

Post-migration:

```
FastAPI → Celery worker → subprocess.Popen("nextflow run X.nf")
                          ↳ Nextflow process runs
                          ↳ Nextflow spawns tool processes
                          ↳ worker tails Nextflow .nextflow.log
                          ↳ parses Nextflow event JSON for progress
```

You add a process layer (the Nextflow runtime) between Celery and the
actual tool. Not slow (10s startup), but more failure modes (Nextflow
runtime crashes, work directory cleanup, `-resume` ambiguity when the
same task is re-submitted by Celery).

### 5. The custom retry semantics differ

Current retry: re-run the whole sample's logfile, check completion
marker. Nextflow retry: per-process `errorStrategy 'retry'` with
`maxRetries`. **Different blast radius**: Nextflow retries one
process; current code retries the entire sample's pipeline. For some
errors (e.g., flaky NFS mounts during alignment) the per-process
retry is better; for others (e.g., upstream tool corrupted the BAM
silently) the whole-sample retry catches it sooner.

### 6. The error pattern detection is broader than `exit_code != 0`

```bash
# Current (preAlignmentQC.sh)
check_logfile "$sample" "fastp" "$log" "fastp_status" "postcheck"
# Internally: greps for tool-specific error patterns AND positive
# completion marker; recovers from cases where exit code is 0 but
# output is corrupt
```

Workflow managers track exit codes only. Migrating loses this defensive
check unless you add explicit `verify:` post-conditions to each
process.

### 7. The `trap_add` composable trap, color logging, and `processbar()`
   helper are all bash idioms that disappear. The web UI doesn't need
   them, but interactive CLI users do.

### 8. `Analysis/` modules are a mix of bash + R + Python

Single-cell scripts under `SingleCellPipe/` use Seurat (R) and
`scanpy`/Scrublet (Python) heavily. Each tool has its own conventions.
Migrating these would mean re-doing each subsystem in
Nextflow/Snakemake's idiom — not a search-and-replace job.

---

## Cost-benefit estimate

| Approach | Effort | Risk | Payoff |
|----------|--------|------|--------|
| **Stay + harden (this PR)** | 2 days | Very low | Real-time progress to UI; safe execution; tractable scripts |
| Containerize per pipeline | 3-5 days | Low | Reproducibility; version pinning; runs identically on dev/prod |
| Migrate `preAlignmentQC` only as PoC | 1 week | Medium | Validate that web integration still works |
| Migrate all to Nextflow | 3-4 weeks | High | DAG/cloud/MultiQC, but lose lab-specific defaults |
| Migrate all to Snakemake | 2-3 weeks | Medium-high | Python idiom matches FastAPI side, but cluster story weaker |

---

## Recommended path forward

Given the FastAPI + Celery web layer is now the user-facing interface,
the most valuable wins are:

1. **Today (this PR)**: Stage A — `set -euo pipefail`, `--dry-run`,
   `::progress::N` protocol + worker integration. Real-time UI updates
   without touching pipeline logic.
2. **Next month**: Stage B — Dockerfile per pipeline, decouple from
   conda. Get reproducibility without Nextflow's runtime overhead.
3. **Later (if needed)**: One pipeline (e.g., `preAlignmentQC`)
   migrated to Nextflow as a parallel option (`/api/v1/pipelines/X?engine=nextflow`).
   Gather data on real-world resume / retry / cluster behavior before
   migrating the rest.
4. **Only if scaling demands**: Full migration. Triggers would be:
   - Need to dispatch jobs to AWS Batch / Google Batch
   - Need to run >10 concurrent pipelines per worker
   - Need cryptographic provenance for clinical reporting
   - Lab grows to a multi-node HPC scheduled environment

Until those triggers fire, the current approach is **simpler, cheaper,
and easier for non-engineers to maintain.** The shell scripts are not
a bug; they are a feature.
