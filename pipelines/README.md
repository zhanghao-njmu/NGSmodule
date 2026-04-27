# pipelines/

This directory holds the new infrastructure for NGSmodule's bash pipelines:

```
pipelines/
├── lib/
│   └── runtime.sh            ← Shared bash helpers (progress, log, dry-run)
└── docker/
    └── preAlignmentQC.Dockerfile   ← Container per pipeline (PoC)
```

The legacy pipeline scripts themselves still live at the repository root
(`GeneralSteps/`, `Analysis/`, `PreparationSteps/`, `PrefetchData/`,
`SingleCellPipe/`) so the existing CLI and conda workflow keep working
unchanged.

## What lives here

### `lib/runtime.sh`

A small bash library scripts can `source` to opt into:

- `set -euo pipefail` (when `NGS_STRICT=1`)
- Structured event protocol (`::progress::N`, `::status::TEXT`,
  `::artifact::kind=X path=Y`, `::metric::K=V`) parsed by the Celery worker
- Coloured CLI logging (`log_info`, `log_warn`, `log_error`, `log_ok`)
- `--dry-run` execution mode (`NGS_DRY_RUN=1`)
- `run_step "label" cmd args...` wrapper combining the above

GeneralSteps/* scripts source this file at the top:

```bash
_ngs_runtime="${shell_folder:-$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/..}/pipelines/lib/runtime.sh"
[[ -f "$_ngs_runtime" ]] && source "$_ngs_runtime"
```

If the file is missing, the script falls back to legacy behaviour — the
hardening is strictly additive.

### `docker/`

Per-pipeline Dockerfiles. Building these and tagging the resulting image
ID into the corresponding `PipelineTemplate` row gives reproducible runs
without rewriting the bash logic.

The `preAlignmentQC.Dockerfile` is the PoC. Replicate the pattern for
`Alignment.Dockerfile`, `postAlignmentQC.Dockerfile`, etc. as you promote
each pipeline.

Build & run:

```bash
docker build \
  -f pipelines/docker/preAlignmentQC.Dockerfile \
  -t ngsmodule/preAlignmentQC:1.0.0 .

docker run --rm \
  -v $(pwd):/opt/ngsmodule:ro \
  -v /data/ngsmodule_work:/data/work \
  -e NGS_TASK_ID=test \
  -e NGS_STRICT=1 \
  ngsmodule/preAlignmentQC:1.0.0 \
  bash /opt/ngsmodule/GeneralSteps/preAlignmentQC.sh ...
```

## Why this layout

Putting hardening helpers in `pipelines/lib/` (instead of mixing with the
existing root-level shell scripts) keeps two things separate:

1. **Lab-specific pipeline logic** — biologists fork, edit, version
2. **Web/runtime contract** — Celery worker depends on the event protocol;
   centralised so a single change updates all pipelines

See [`../docs/PIPELINE_MIGRATION_ANALYSIS.md`](../docs/PIPELINE_MIGRATION_ANALYSIS.md)
for the deeper architectural reasoning, including the explicit decision
*not* to migrate to Snakemake/Nextflow at this stage.
