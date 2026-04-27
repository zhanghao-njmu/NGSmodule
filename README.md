# 🧬 NGSmodule

**Self-hosted next-generation sequencing analysis platform — from FASTQ to insight.**

<div align="center">

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/fastapi-0.109-009688.svg)](https://fastapi.tiangolo.com/)
[![React](https://img.shields.io/badge/react-18-61dafb.svg)](https://react.dev/)
[![PostgreSQL](https://img.shields.io/badge/postgres-15-336791.svg)](https://www.postgresql.org/)
[![Docker](https://img.shields.io/badge/docker-ready-2496ed.svg)](https://www.docker.com/)

[Quickstart](#-quickstart) ·
[Pipelines](#-supported-pipelines) ·
[Architecture](#-architecture) ·
[For wet-lab biologists](#-for-wet-lab-biologists) ·
[For bioinformaticians](#-for-bioinformaticians) ·
[Deployment](#-deployment)

</div>

---

## What is NGSmodule?

NGSmodule turns your sequencing data — **FASTQ, BAM, VCF** — into reproducible
analyses without making researchers write shell scripts. It runs the standard
NGS toolchain (`fastp`, `BWA-MEM2`, `STAR`, `HISAT2`, `Salmon`, `Kallisto`,
`samtools`, `GATK`, …) under the hood, but exposes:

- A **web UI** for project/sample/file management with batch CSV import
- A **template-driven pipeline runner** with parameter recommendation
- A **real-time task monitor** (WebSocket) so you can watch a 12-hour
  alignment from your phone
- An **AI assistant** (Claude/OpenAI/local) that reads your QC reports and
  suggests parameter changes
- **Multi-user RBAC** with per-user storage quotas — suitable for shared
  lab compute servers

If you've ever maintained a Snakemake/Nextflow pipeline for non-coders,
NGSmodule is the missing self-service layer on top.

---

## 🧪 Supported Pipelines

Built-in templates target common NGS workflows. Each is a versioned
`PipelineTemplate` you can clone, edit parameters on, and execute against
any subset of samples.

| Category | Template | Tools | Typical use |
|----------|----------|-------|-------------|
| **Quality Control** | Pre-Alignment QC | `fastp`, `FastQC`, `Trimmomatic` | Adapter trimming, quality filtering before mapping |
| | Post-Alignment QC | `RSeQC`, `Qualimap`, `dupRadar` | Strand specificity, duplication, coverage |
| **Alignment** | Sequence Alignment | `BWA-MEM2`, `STAR`, `HISAT2`, `samtools` | DNA-seq / RNA-seq mapping with stats |
| **RNA-seq** | Gene Expression Quantification | `Salmon`, `Kallisto`, `featureCounts` | TPM / counts matrices |
| | Differential Expression | `DESeq2`, `edgeR`, `limma-voom` | Two-condition comparisons + volcano plot |
| **DNA-seq** | GATK Germline Variant Calling | `GATK4 HaplotypeCaller`, `bcftools` | Single-sample / joint variant calling |
| | GATK Somatic Variant Calling | `Mutect2`, `GATK4` | Tumor-normal somatic mutations |
| | Copy Number Variation | `CNVkit` | Targeted/WES CNV |

> Pipelines are shell scripts under `GeneralSteps/` — fully readable, fully
> hackable. You can register new templates without touching the backend code.

### Sample → Result data flow

```
FASTQ (uploaded)
    │
    ▼
┌───────────────────┐    ┌───────────────────┐    ┌───────────────────┐
│ Pre-Alignment QC  │ →  │   Alignment       │ →  │ Post-Alignment QC │
│ (fastp / FastQC)  │    │ (STAR / BWA-MEM2) │    │ (Qualimap, dupRadar)│
└───────────────────┘    └───────────────────┘    └───────────────────┘
                                   │
                                   ├─→ RNA-seq quantification (Salmon/Kallisto)
                                   │      └─→ Differential expression (DESeq2)
                                   │
                                   └─→ Variant calling (GATK)
                                          └─→ CNV analysis (CNVkit)
```

Each box becomes a `Task` in the database with its own log file, progress
bar, and result artefacts (FASTQ, BAM, VCF, count matrices, HTML reports).

---

## 🚀 Quickstart

### Single-host development

```bash
git clone https://github.com/zhanghao-njmu/NGSmodule.git
cd NGSmodule

# Bring up postgres + redis + minio + backend + worker + frontend
docker compose up -d

# Initialise database schema and seed the 8 built-in pipeline templates
docker compose exec backend alembic upgrade head
docker compose exec backend python init_pipeline_templates.py
docker compose exec backend python create_admin.py  # creates admin/admin123
```

Open <http://localhost:3000>.
API docs at <http://localhost:8000/api/v1/docs>.

### First analysis (5 minutes)

1. **Login** as `admin` (change the password in *Settings* immediately)
2. *Projects* → *New Project* → "RNA-seq pilot study"
3. *Samples* → *Import CSV* (use the template at
   [`docs/sample-template.csv`](docs/sample-template.csv); columns:
   `sample_id, run_id, group_name, layout, batch_id`)
4. *Files* → *Upload FASTQ* (drag-and-drop, supports up to 50 GB per file
   via MinIO multipart)
5. *Pipelines* → *Pre-Alignment QC* → *Execute* → select samples
6. *Tasks* → watch live progress (WebSocket, no refresh needed)
7. *Results* → interactive Plotly QC report (per-base quality, GC content,
   duplication, adapter content)

---

## 🏗️ Architecture

```
                 ┌────────────────────────────────────┐
                 │  Browser (SPA, rsbuild + React 18) │
                 │  TanStack Query · AntD · Plotly    │
                 └─────────┬──────────────┬───────────┘
                           │ HTTPS        │ WSS
                           ▼              ▼
                 ┌────────────────────────────────────┐
                 │      Nginx (TLS, static, proxy)    │
                 └─────────┬──────────────┬───────────┘
                           ▼              ▼
                 ┌────────────────────────────────────┐
                 │   FastAPI (uvicorn / ASGI)         │
                 │   ~146 REST + 1 WebSocket          │
                 │   JWT auth · slowapi rate-limit    │
                 └──┬────┬────┬────┬────┬─────────────┘
                    │    │    │    │    │
        ┌───────────┘    │    │    │    └──────────────┐
        ▼                ▼    ▼    ▼                   ▼
   ┌─────────┐     ┌────────┐ │ ┌────────┐        ┌─────────┐
   │ Postgres│     │ Redis  │ │ │ MinIO  │        │ Sentry  │
   │ (ORM    │     │ broker+│ │ │ S3 obj │        │ +       │
   │  state) │     │ pubsub │ │ │ store  │        │ Prom.   │
   └─────────┘     └────┬───┘ │ │ FASTQ/ │        │ (opt.)  │
                        ▼     │ │ BAM/VCF│        └─────────┘
                  ┌─────────────┐ └────────┘
                  │ Celery      │
                  │ workers ×N  │  (4 queues)
                  │  default    │   - notifications, audit
                  │  admin      │   - cleanup, sync
                  │  backup     │   - pg_dump, file tar
                  │  ngs_pipeline   - actual NGS jobs
                  └─────────────┘
```

**Key design decisions**

| Component | Choice | Why |
|-----------|--------|-----|
| Backend framework | FastAPI | OpenAPI-native, async, Pydantic validation |
| Workers | Celery + Redis | Multi-queue split (heavy NGS jobs vs. light admin tasks); horizontally scalable |
| Object storage | MinIO (S3 API) | Self-hosted, handles GB-scale FASTQ uploads via multipart presigned URLs |
| Real-time | Redis pub/sub → WebSocket | Cross-replica fanout; auto-reconnect with exponential backoff |
| Frontend bundler | rsbuild (Rspack) | Rust-based; 5–10× faster than webpack/vite for large SPAs |
| Server-state | TanStack Query | 74 typed hooks, optimistic mutations, realtime cache invalidation |
| Charts | Plotly only | Genomics-grade visualizations (volcano plots, Manhattan plots, heatmaps) |
| AI provider | Pluggable (mock / Claude / OpenAI) | `AI_PROVIDER` env var hot-swaps backend; falls back to mock when offline |
| Observability | Sentry + Prometheus + audit log | Optional; gracefully disabled when libs missing |

See [docs/MULTI_REPLICA_DEPLOYMENT.md](docs/MULTI_REPLICA_DEPLOYMENT.md) for
horizontal scaling and Kubernetes manifests.

---

## 🧫 For wet-lab biologists

You shouldn't need to read a single line of code:

- **Drag-and-drop FASTQ uploads** — server resumes interrupted transfers
- **CSV sample sheets** — paste from Excel, system validates and reports
  duplicates / missing fields
- **Run any built-in pipeline** by clicking *Execute* and picking samples
- **AI parameter recommendations** — click "Get Recommendations" in the
  pipeline modal; the assistant analyses your past successful runs (or
  consults Claude/OpenAI if configured) and pre-fills sensible parameters
  with a confidence score
- **Browser-based result viewer** — every QC report renders as interactive
  Plotly charts; volcano plots are zoomable, gene names searchable
- **Real-time progress** — start a 12-hour alignment, close your laptop,
  re-open it tomorrow → progress bar is exactly where it should be
- **Notifications** — email + in-app when tasks finish/fail
- **AI assistant** — ask questions like "why does sample S07 have 18%
  duplication?" in the chat; it reads your QC results and answers in plain
  language

---

## 🧑‍💻 For bioinformaticians

If you live on the command line, NGSmodule still respects you:

- **Every pipeline is a shell script** in `GeneralSteps/` — readable, forkable,
  testable
- **Pipeline templates are JSON-defined** parameter schemas in
  `init_pipeline_templates.py`; add new ones without touching Python
- **Full REST API** at `/api/v1/docs` — automate the entire UI, including
  CSV import, batch execution, result download
- **Background tasks via Celery** — drop a `@celery_app.task` into
  `app/workers/`, route to a queue (default / admin / backup / ngs_pipeline)
- **Database migrations** managed by Alembic
- **Prometheus `/metrics`** for FastAPI request timings, plus per-task Celery
  metrics through Flower (`:5555`)
- **Audit log** records every admin action (user role change, password
  reset, system config update) with IP and User-Agent — exportable as
  JSON/CSV for compliance audits
- **Self-hostable AI** — point `AI_PROVIDER=claude` + `ANTHROPIC_API_KEY` at
  Anthropic's API (with **prompt caching** enabled to keep costs low), or
  swap in a local model by implementing the `AIProvider` interface in
  `app/services/ai_providers/`

### REST API examples

```bash
# Login
TOKEN=$(curl -s -X POST http://localhost:8000/api/v1/auth/login \
  -H "Content-Type: application/json" \
  -d '{"username":"admin","password":"admin123"}' | jq -r .access_token)

# Batch import 96-well sample sheet
curl -X POST http://localhost:8000/api/v1/samples/import \
  -H "Authorization: Bearer $TOKEN" \
  -F "project_id=$PROJECT_ID" \
  -F "file=@plate_001.csv"

# Execute Pre-Alignment QC on every sample in a project
curl -X POST http://localhost:8000/api/v1/pipelines/batch-execute \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "template_id": "preAlignmentQC",
    "project_id": "'$PROJECT_ID'",
    "sample_ids": ["s1", "s2", "..."],
    "task_name_prefix": "QC pilot",
    "parameters": {"min_quality": 20, "min_length": 50}
  }'
```

---

## 🚢 Deployment

### Single-node (Docker Compose)

```bash
docker compose -f docker-compose.prod.yml up -d
```

Suitable for ≤ 10 concurrent users / labs.

### Multi-replica (production)

```bash
# Light queue × 4, heavy NGS queue × 2, beat scheduler × 1
docker compose \
  -f docker-compose.prod.yml \
  -f docker-compose.scale.yml \
  up -d \
  --scale celery-worker=4 \
  --scale celery-worker-heavy=2
```

Recommended hardware sizing:

| Workload | Replicas (default) | Replicas (heavy) | Resource per heavy worker |
|----------|--------------------|------------------|---------------------------|
| Single lab (1–5 users) | 1 | 1 | 8 vCPU, 16 GB RAM |
| Department (5–25 users) | 4 | 2 | 16 vCPU, 64 GB RAM |
| Core facility (50+ users) | 8 | 4 | 32 vCPU, 256 GB RAM |

### Kubernetes

See [docs/MULTI_REPLICA_DEPLOYMENT.md](docs/MULTI_REPLICA_DEPLOYMENT.md) for
manifests covering `Deployment`, `Service`, `HorizontalPodAutoscaler`, and
`PersistentVolumeClaim` for FASTQ storage.

### Optional integrations

| Feature | How to enable |
|---------|---------------|
| Sentry error tracking | `SENTRY_DSN=https://...sentry.io/...` |
| Prometheus metrics | `pip install prometheus-fastapi-instrumentator` (auto-mounted at `/metrics`) |
| Real Claude AI | `AI_PROVIDER=claude` + `ANTHROPIC_API_KEY=...` |
| Real OpenAI | `AI_PROVIDER=openai` + `OPENAI_API_KEY=...` |
| Multi-replica realtime | `REDIS_URL=redis://redis:6379/0` (defaults to in-memory) |
| Backup retention | `BACKUP_DIR=/data/backups`, `BACKUP_RETENTION_DAYS=30` |

---

## 🔒 Security

- **Authentication**: JWT (HS256) with refresh tokens; bcrypt password hashing
- **Authorization**: User / Admin roles + per-resource ownership checks
- **Rate limiting**: per-IP and per-user, Redis-backed (auto-disabled in tests)
- **Audit log**: every privileged action recorded with IP + User-Agent + Request ID
- **Data isolation**: hard tenancy by `user_id` on every read query
- **Soft delete**: deletion of projects/samples is reversible by admins for 30 days
- **Backups**: scheduled `pg_dump` + tar of file storage, SHA-256 checksum,
  configurable retention

---

## 📊 Status

| Layer | Coverage |
|-------|----------|
| Backend API | 146 routes across 14 modules |
| Backend tests | 33/33 passing |
| Frontend pages | 18 pages, all on TanStack Query |
| Frontend type-check | 0 errors |
| Frontend lint | 0 errors |
| Production build | rsbuild · 1.9 MB gzipped · 8s |
| AI Intelligence | 28 endpoints (mock by default; Claude/OpenAI ready) |
| Admin Enhanced | 27 endpoints (audit, alerts, backups, jobs, metrics) |

---

## 📚 Documentation

- [Multi-Replica Deployment Guide](docs/MULTI_REPLICA_DEPLOYMENT.md)
- [Frontend-Backend Integration Report](FRONTEND_BACKEND_INTEGRATION_REPORT.md)
- [Production Deployment](PRODUCTION_DEPLOYMENT.md)
- [Security Best Practices](SECURITY_BEST_PRACTICES.md)
- [API Test Collection](API_TEST_COLLECTION.md)
- [Monitoring & Logging Guide](MONITORING_LOGGING_GUIDE.md)
- Live API docs: <http://localhost:8000/api/v1/docs>

---

## 🤝 Contributing

NGSmodule is built primarily as a self-hosted tool, but contributions are
welcome. Common areas:

- **New pipeline templates** — drop a shell script into `GeneralSteps/`,
  register it in `init_pipeline_templates.py`
- **Additional AI providers** — implement `AIProvider` from
  `app/services/ai_providers/base.py`
- **Bioinformatics-specific visualizations** — Plotly volcano/Manhattan/heatmap
  variants live under `frontend/src/components/charts/`

Run the test suite before opening a PR:

```bash
# Backend
cd backend && pytest

# Frontend
cd frontend && npm run validate   # tsc + eslint + prettier
```

---

## 📄 License

MIT — see [LICENSE](LICENSE).

If NGSmodule is useful in your research, please cite this repository in your
methods section.
