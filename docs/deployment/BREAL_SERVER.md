# Deploying NGSmodule on the BREAL Server

This is the host-specific recipe for `BREAL-server` (zhanghao). It
encodes the conventions already used by all the other services on this
machine:

- All docker services live under `/ssd/docker/services/<name>/` with a
  thin `docker-compose.yml` + bind-mounted `data/` directory.
- nginx-proxy-manager (in `/ssd/docker/services/nginx-proxy-manager/`)
  is the SSL/HTTPS terminator — services bind 127.0.0.1-only ports and
  let NPM proxy `:80/:443 → 127.0.0.1:NNNN`.
- NGS data follows a fixed layout:
  - `/data/lab/{user}/{datatype}/{project}/rawdata/` — raw fastq.gz
    (cold, on the 29T spinning disk)
  - `/ssd/reference/` — shared iGenomes / SortmeRNA / FastQ_Screen
  - `/ssd/lab/{user}/{datatype}/{project}/` — hot working tree
    (NGSmodule_work + NGSmodule_analysis + multiqc + Sample_info.csv)

## Port budget

| Range | Used by |
|---|---|
| 3000–3003 | openwebui / docker-nginx |
| 4000–4005 | portainer / uptime-kuma / tika / libretranslate |
| 5003–5212 | dify-plugin-daemon / cloudreve |
| 8000–8082 | jupyterhub / searxng / grobid / paddleocr |
| 9099 / 9898 | openwebui-pipelines / backrest |
| 11434 | ollama |

**Allocation for NGSmodule** (chosen to avoid collisions, all
localhost-only — NPM publishes to public via 80/443):

| Port | Service |
|---|---|
| 4173 | frontend (nginx) |
| 4180 | backend FastAPI |
| 4181 | flower (Celery monitoring, optional) |

Internal services (postgres / redis / minio) do **not** publish a host
port — they're reachable only inside the docker network.

## One-time setup

```bash
# 1. Service home (sibling to openwebui/, dify/, …)
sudo install -d -o "$USER" -g "$USER" /ssd/docker/services/ngsmodule
cd /ssd/docker/services/ngsmodule

# 2. Persistent data + backup directories
mkdir -p data/{postgres,redis,minio,uploads,storage} backups/postgres logs

# 3. Bind config from the source repo (so `git pull` updates compose
#    files in place; no copy/diff to drift).
ln -sf /home/zhanghao/programs/NGSmodule/docker-compose.prod.yml ./docker-compose.prod.yml
ln -sf /home/zhanghao/programs/NGSmodule/docker-compose.host.yml ./docker-compose.host.yml

# 4. Generate strong secrets and write .env.
cat > .env <<EOF
# === MUST-SET ===
SECRET_KEY=$(openssl rand -hex 32)
JWT_SECRET_KEY=$(openssl rand -hex 32)
POSTGRES_USER=ngsmodule
POSTGRES_PASSWORD=$(openssl rand -hex 24)
POSTGRES_DB=ngsmodule
REDIS_PASSWORD=$(openssl rand -hex 24)
MINIO_ROOT_USER=$(openssl rand -hex 8)
MINIO_ROOT_PASSWORD=$(openssl rand -hex 24)
MINIO_ACCESS_KEY=\${MINIO_ROOT_USER}
MINIO_SECRET_KEY=\${MINIO_ROOT_PASSWORD}
FLOWER_PASSWORD=$(openssl rand -hex 16)

# === Image tag ===
IMAGE_TAG=latest

# === NGS pipeline ===
NGS_PIPELINE_DIR=/app/NGSmodule
LCBIO_HOST_DIR=/opt/lcbio    # adjust if you installed lcbio elsewhere

# === CORS / public hostname (set after NPM domain is wired) ===
BACKEND_CORS_ORIGINS=http://localhost:4173

# === Frontend bake-in ===
VITE_API_URL=/api/v1
VITE_WS_URL=/api/v1/ws
EOF

chmod 600 .env

# 5. First start
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml \
  --env-file .env up -d --build

# 6. Wait for the database, run migrations
until docker inspect ngsmodule-postgres --format '{{.State.Health.Status}}' | grep -q healthy; do sleep 3; done
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend \
  alembic upgrade head

# 7. Verify
curl -s http://127.0.0.1:4180/health           # {"status":"healthy"}
curl -sI http://127.0.0.1:4173/ | head -1      # HTTP 200
```

## Bootstrap admin account

```bash
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend python -c "
from app.core.database import SessionLocal
from app.models.user import User
from app.core.security import get_password_hash
import uuid
db = SessionLocal()
u = User(
    id=uuid.uuid4(),
    username='admin',
    email='admin@breal.local',
    hashed_password=get_password_hash('CHANGE_ME_FIRST_LOGIN'),
    role='admin',
    is_active=True,
    storage_quota=10**12,
)
db.add(u); db.commit()
print('Admin created:', u.id)
"
```

Login at `http://<server-ip>:4173` (or your NPM domain), change the
admin password from the UI immediately.

## Wire HTTPS via nginx-proxy-manager

```bash
# 1. Bring NPM up (configured but not running by default).
cd /ssd/docker/services/nginx-proxy-manager
docker compose up -d

# 2. Visit http://<server-ip>:81
#    Default login: admin@example.com / changeme  → change it.

# 3. Add a Proxy Host:
#      Domain Names:    ngs.your-domain.com
#      Forward Hostname: 127.0.0.1
#      Forward Port:     4173
#      Cache Assets, Block Common Exploits, Websockets Support — all on.
#    SSL tab:
#      Request a new SSL Certificate — Let's Encrypt → Force SSL.
```

After the cert is issued, update `.env`:

```bash
sed -i 's|^BACKEND_CORS_ORIGINS=.*|BACKEND_CORS_ORIGINS=https://ngs.your-domain.com|' .env
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml \
  --env-file .env up -d --no-deps backend
```

## Verifying NGS data wiring

The host bind mounts let the backend reach the canonical layout
`/data/lab/{user}/...` (raw) and `/ssd/lab/{user}/...` (work). To
register an existing project for a user, the simplest sanity check is:

```bash
# Inside the backend container, the same paths the user sees on the
# host should be visible read/write.
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend ls /ssd/lab/caiwenjin/RNA-seq/aging_msc/NGSmodule_work | head
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend ls /data/lab/caiwenjin/RNA-seq/aging_msc/rawdata | head
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend ls /ssd/reference/iGenomes | head
```

If those commands print the user's actual files, the path-mapping is
correct and the API can scan them via the standard NGSmodule pipeline
config.

## Common operations

```bash
# Stop the stack (keep data)
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml down

# Rolling update after `git pull`
cd /home/zhanghao/programs/NGSmodule && git pull
cd /ssd/docker/services/ngsmodule
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml \
  --env-file .env build backend celery-worker frontend
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml \
  --env-file .env up -d
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec backend \
  alembic upgrade head

# Postgres backup (daily cron candidate)
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml exec -T postgres \
  pg_dump -U ngsmodule ngsmodule | gzip > backups/postgres/$(date +%F).sql.gz

# Tail backend / worker logs
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml logs -f backend
docker compose -f docker-compose.prod.yml -f docker-compose.host.yml logs -f celery-worker

# Disk usage by NGSmodule volumes
du -sh /ssd/docker/services/ngsmodule/data/* 2>/dev/null
```

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `port already in use` on first boot | Another service grabbed 4173/4180/4181 | `ss -ltn | grep -E ':41(7[3-9]|8[01])'`, pick a different port in `docker-compose.host.yml`. |
| Backend can't see `/ssd/lab/...` | `:rw` bind mount didn't propagate | Restart with `--force-recreate backend`. |
| `Permission denied` writing to data | docker user mismatch on /ssd/lab | The mounts inherit host UID/GID; run as your normal user, not root. The work dir at `/ssd/lab/<user>` should already be group-writable by `shared`. |
| `failed to bind host port` | Stale port-mapping cache | `docker compose ... down && up -d --force-recreate`. |
| `alembic upgrade` errors | Image rebuilt without DB ready | Wait until `docker inspect ngsmodule-postgres ... healthy` before running migrations. |
