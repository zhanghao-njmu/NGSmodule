# Multi-Replica Celery Deployment Guide

This document describes how to scale the NGSmodule Celery workers
horizontally for production workloads.

## Architecture

NGSmodule uses **four task queues**:

| Queue          | Purpose                                       | Latency        |
|----------------|-----------------------------------------------|----------------|
| `default`      | Notifications, audit log writes               | < 100ms        |
| `admin`        | System health checks, cleanup, job sync       | < 10s          |
| `backup`       | Database / file backups                       | minutes – hours|
| `ngs_pipeline` | Sequence alignment, quantification, etc.      | hours – days   |

Routing is configured in `app/workers/celery_app.py` via `task_routes`.

## Worker Pools

Two distinct worker pools should be deployed:

### 1. Default workers (light, scaled horizontally)

Process `default` and `admin` queues. Horizontally scalable; each replica
should be small (≤ 4 vCPU, ≤ 2GB).

```bash
celery -A app.workers.celery_app worker \
  --concurrency=4 \
  --max-tasks-per-child=100 \
  --queues=default,admin \
  --hostname=worker-default@%h
```

### 2. Heavy workers (NGS + backup, vertically scaled)

Process `ngs_pipeline` and `backup` queues. Few replicas, each with
generous resources.

```bash
celery -A app.workers.celery_app worker \
  --concurrency=2 \
  --max-tasks-per-child=10 \
  --queues=ngs_pipeline,backup \
  --hostname=worker-heavy@%h \
  --time-limit=86400 \
  --soft-time-limit=82800
```

### 3. Beat scheduler (single replica only)

The Celery beat scheduler must run as a single replica to avoid
duplicate periodic task execution:

```bash
celery -A app.workers.celery_app beat --loglevel=info
```

## Docker Compose Deployment

Use the provided `docker-compose.scale.yml` overlay:

```bash
# Set scaling factors
export CELERY_DEFAULT_REPLICAS=4
export CELERY_HEAVY_REPLICAS=2
export CELERY_DEFAULT_CONCURRENCY=4
export CELERY_HEAVY_CONCURRENCY=2

# Start with multi-replica config
docker compose \
  -f docker-compose.prod.yml \
  -f docker-compose.scale.yml \
  up -d
```

Or scale at runtime:

```bash
docker compose -f docker-compose.prod.yml up -d \
  --scale celery-worker=4
```

## Docker Swarm Deployment

```bash
docker stack deploy \
  -c docker-compose.prod.yml \
  -c docker-compose.scale.yml \
  ngsmodule
```

## Kubernetes Deployment

For Kubernetes, deploy each worker pool as a separate `Deployment`:

```yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: celery-worker-default
spec:
  replicas: 4
  selector:
    matchLabels: { app: celery-worker-default }
  template:
    metadata:
      labels: { app: celery-worker-default }
    spec:
      containers:
        - name: celery
          image: ngsmodule/backend:latest
          command: ["celery", "-A", "app.workers.celery_app", "worker"]
          args:
            - "--concurrency=4"
            - "--queues=default,admin"
            - "--max-tasks-per-child=100"
          resources:
            requests:
              cpu: "500m"
              memory: "512Mi"
            limits:
              cpu: "2"
              memory: "2Gi"
---
apiVersion: apps/v1
kind: Deployment
metadata:
  name: celery-worker-heavy
spec:
  replicas: 2
  selector:
    matchLabels: { app: celery-worker-heavy }
  template:
    metadata:
      labels: { app: celery-worker-heavy }
    spec:
      containers:
        - name: celery
          image: ngsmodule/backend:latest
          command: ["celery", "-A", "app.workers.celery_app", "worker"]
          args:
            - "--concurrency=2"
            - "--queues=ngs_pipeline,backup"
            - "--max-tasks-per-child=10"
          resources:
            requests:
              cpu: "4"
              memory: "8Gi"
            limits:
              cpu: "8"
              memory: "16Gi"
---
# beat must be exactly 1 replica
apiVersion: apps/v1
kind: Deployment
metadata:
  name: celery-beat
spec:
  replicas: 1
  strategy: { type: Recreate }  # never run two beats simultaneously
  selector:
    matchLabels: { app: celery-beat }
  template:
    metadata:
      labels: { app: celery-beat }
    spec:
      containers:
        - name: celery
          image: ngsmodule/backend:latest
          command: ["celery", "-A", "app.workers.celery_app", "beat"]
```

## Worker Hardening (already enabled)

The following safety settings are configured in `celery_app.py` for
multi-replica deployments:

- `worker_prefetch_multiplier=1` — prevents a single worker from hoarding
  long-running tasks
- `task_acks_late=True` — tasks are only acknowledged after successful
  completion, so a crashed worker's tasks get redelivered
- `task_reject_on_worker_lost=True` — same as above for SIGKILL/OOM
- `broker_connection_retry_on_startup=True` — graceful Redis reconnection

## Monitoring

- **Flower** (already in compose): real-time dashboard at `:5555`
- **Prometheus** (P1-2): worker metrics auto-exposed at `/metrics`
- **Sentry** (P1-3): worker exceptions automatically captured via the
  CeleryIntegration

## Scaling Heuristics

| Workload                  | Default replicas | Heavy replicas |
|---------------------------|------------------|----------------|
| Dev / single-tenant       | 1                | 1              |
| 1 – 10 active users       | 2                | 1              |
| 10 – 50 active users      | 4                | 2              |
| 50+ active users          | 8                | 4              |
| Burst (batch processing)  | 8 – 16           | 4 – 8          |
