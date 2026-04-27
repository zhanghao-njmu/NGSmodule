# NGSmodule documentation

This directory holds all long-form documentation. The repository root
[`README.md`](../README.md) is the entry point; everything else lives
here, organised by audience.

## Map

```
docs/
├── deployment/         ← Production operators
│   ├── GUIDE.md            – step-by-step deployment
│   ├── CHECKLIST.md        – pre-flight checklist
│   ├── PRODUCTION.md       – production-grade hardening
│   ├── MONITORING.md       – Sentry / Prometheus / logging setup
│   ├── MULTI_REPLICA.md    – multi-replica Celery deployment
│   └── SECURITY.md         – security best practices
│
├── development/        ← Contributors
│   ├── MASTER_PLAN.md
│   ├── COMPREHENSIVE_DEVELOPMENT_PLAN.md
│   ├── PLAN.md
│   ├── IMPLEMENTATION_ROADMAP.md
│   ├── EXECUTIVE_SUMMARY.md
│   ├── BACKEND_API_IMPLEMENTATION_PLAN.md
│   ├── BACKEND_OPTIMIZATION_PLAN.md
│   └── FRONTEND_BACKEND_INTEGRATION_REPORT.md
│
├── testing/            ← Test authors
│   ├── BACKEND_API_TESTING_GUIDE.md
│   ├── E2E_GUIDE.md
│   ├── FRONTEND_CHECKLIST.md
│   └── INTEGRATION.md
│
├── api/                ← API consumers
│   └── TEST_COLLECTION.md  – curl / Postman examples
│
├── PIPELINE_MIGRATION_ANALYSIS.md   ← Snakemake/Nextflow trade-off study
└── sample-template.csv               ← Minimal sample sheet for the UI
```

## Where to start

| If you want to… | Read |
|-----------------|------|
| Deploy NGSmodule on a server | [`deployment/GUIDE.md`](deployment/GUIDE.md) |
| Scale to multiple replicas | [`deployment/MULTI_REPLICA.md`](deployment/MULTI_REPLICA.md) |
| Harden for production | [`deployment/PRODUCTION.md`](deployment/PRODUCTION.md) + [`deployment/SECURITY.md`](deployment/SECURITY.md) |
| Set up monitoring | [`deployment/MONITORING.md`](deployment/MONITORING.md) |
| Add a new backend endpoint | [`development/BACKEND_API_IMPLEMENTATION_PLAN.md`](development/BACKEND_API_IMPLEMENTATION_PLAN.md) |
| Write a frontend test | [`testing/FRONTEND_CHECKLIST.md`](testing/FRONTEND_CHECKLIST.md) |
| Run end-to-end tests | [`testing/E2E_GUIDE.md`](testing/E2E_GUIDE.md) |
| Hit the API directly | [`api/TEST_COLLECTION.md`](api/TEST_COLLECTION.md) |
| Decide whether to migrate pipelines | [`PIPELINE_MIGRATION_ANALYSIS.md`](PIPELINE_MIGRATION_ANALYSIS.md) |

## Living docs

- API reference (auto-generated from FastAPI): <http://localhost:8000/api/v1/docs>
- AI providers reference: see `backend/app/services/ai_providers/base.py`
