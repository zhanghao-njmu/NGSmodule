# pipelines/

The new pipeline framework for NGSmodule. Sample-first, dependency-resolving,
schema-validated, cluster-aware, dry-runnable, scope:project-aware.

```
pipelines/
├── core/                  ← One directory per pipeline (18 today)
│   └── <Name>/
│       ├── meta.yml       schema, requires, resources, scope
│       ├── env.yml        conda environment
│       ├── pipeline.sh    run_<Name>_for_{sample,project}() function
│       ├── scripts/       (optional) bundled R / python helpers
│       └── tests/         dry-run fixture
├── lib/                   ← framework primitives (orchestrator, schema, ...)
├── modules/               ← reusable tool wrappers (18 today)
├── profiles/              ← cluster scheduler adapters (slurm/sge/lsf/local)
└── docker/                ← per-pipeline Dockerfiles (PoC)
```

## Quickstart

```bash
ngsmodule list                                        # all available pipelines
ngsmodule deps DifferentialExpression                 # see what 'run' would resolve
ngsmodule run Alignment -c project.cfg                # execute (auto-resolves prereqs)
ngsmodule resume Alignment -c project.cfg             # re-run only failed samples
ngsmodule status -c project.cfg                       # sample × stage grid
ngsmodule report -c project.cfg -o run.html           # self-contained HTML report
ngsmodule init MyPipe --requires Alignment            # scaffold a new pipeline
ngsmodule lint                                        # offline schema validation
ngsmodule test                                        # dry-run regression suite
```

## Pipeline catalogue (18)

| Scope     | Pipeline                  | Category        | Tool(s)                                         |
| --------- | ------------------------- | --------------- | ----------------------------------------------- |
| sample    | preAlignmentQC            | QC              | fastp, fastqc                                   |
| sample    | Alignment                 | Alignment       | STAR / bwa-mem2 / hisat2 / bismark              |
| sample    | postAlignmentQC           | QC              | RSeQC, preseq, mosdepth, goleft                 |
| sample    | Quantification            | Quantification  | featureCounts                                   |
| sample    | MethylationExtraction     | Methylation     | bismark methylation_extractor                   |
| sample    | CircRNA                   | Quantification  | STAR + CIRCexplorer2                            |
| sample    | VariantCalling            | VariantCalling  | bcftools mpileup/call                           |
| sample    | GATK_germline_short       | VariantCalling  | GATK4 HaplotypeCaller                           |
| sample    | GATK_somatic_short        | VariantCalling  | GATK4 Mutect2                                   |
| sample    | GATK_CNV                  | VariantCalling  | GATK4 ModelSegments (consumes PoN)              |
| sample    | Strelka2_germline         | VariantCalling  | Strelka2                                        |
| sample    | Strelka2_somatic          | VariantCalling  | Strelka2 T/N                                    |
| project   | MergeCounts               | Quantification  | awk merge per-sample counts                     |
| project   | postQuantificationQC      | Quantification  | base-R PCA + correlation                        |
| project   | BatchCorrection           | Quantification  | sva ComBat-seq / limma                          |
| project   | DifferentialExpression    | Quantification  | edgeR (auto-uses corrected matrix)              |
| project   | DifferentialMethylation   | Methylation     | base-R CpG Welch t-test                         |
| project   | GATK_CNV_PoN              | VariantCalling  | GATK4 CreateReadCountPanelOfNormals             |

## What you get for free

- **State machine.** `<work_dir>/<sample>/.state.json` records each stage's
  status, started_at, completed_at, params, tools, inputs, outputs.
- **Auto-discovery + dependency resolution.** Declaring `requires:` in
  meta.yml is enough — the orchestrator transitively resolves prereqs.
- **Project-scope.** Set `scope: project` and define
  `run_<Name>_for_project()`. State goes to `_project/.state.json`.
- **Schema validation.** Typed `params_schema:` (enum, int, float,
  bool, path, string with pattern) is checked before any sample runs.
- **Cluster execution.** `--profile slurm|sge|lsf|local` swaps the
  per-sample submit logic without touching pipeline code.
- **Dry-run + tests + CI.** Every pipeline ships a `tests/` fixture;
  `ngsmodule test` runs the lot in dry-run mode (no real tools needed).
- **HTML report.** Self-contained, no JS, no CDN. Status grid + Gantt
  + resource peaks + tool tally + failure breakdown + lineage.

## Migration & contribution

For a step-by-step guide to porting a legacy script into a new pipeline,
see [`../docs/MIGRATION_GUIDE.md`](../docs/MIGRATION_GUIDE.md). The
fastest start is `ngsmodule init <Name>`.

## Out of scope

`SingleCellPipe/` (cellranger / Seurat / SCP custom R package) is a
separate ~11k-line subframework and has not been migrated. See the
"Migration coverage" section of MIGRATION_GUIDE.md for an outline of
what porting it would involve.
