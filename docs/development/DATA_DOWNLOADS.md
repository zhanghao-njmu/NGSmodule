# Vendor Data Downloads

NGSmodule has a unified API for fetching sequencing deliveries from
external partners (联川生物 / 诺禾致源 / ...) directly into the platform.
The frontend talks to this API; users no longer need to ssh in and run
the vendor's CLI manually.

## Architecture

```
Frontend ──HTTP──▶ Backend (FastAPI)
                       │
                       ├── data_download_service ──▶ DataProviderBase
                       │                              ├── LcBioProvider  (subprocess wraps lcbio bash CLI)
                       │                              └── NovogeneProvider (stub)
                       │
                       └── Celery worker (data_downloads.watch_progress)
                                  │
                                  ├── tails vendor log file
                                  ├── updates DownloadJob row
                                  └── publishes realtime:task:<id> events
                                                │
Frontend ◀──WebSocket────────────────────────────
```

## Endpoints

All under `/api/v1/data-downloads/` unless noted; require a logged-in user.

| Method | Path | Purpose |
|--------|------|---------|
| GET    | `/vendors`                     | list supported vendor ids |
| GET    | `/sessions/{vendor}`           | is the vendor daemon up? |
| POST   | `/sessions`                    | open a session (email+pw OR `credential_id`) |
| DELETE | `/sessions/{vendor}`           | kill the vendor daemon |
| POST   | `/jobs`                        | start a download |
| GET    | `/jobs`                        | list current user's jobs (refreshes in-flight) |
| GET    | `/jobs/{id}`                   | single job (re-parses log) |
| DELETE | `/jobs/{id}`                   | cancel (lcbio: kills daemon, affects all in-flight) |

Saved-credentials API (under `/api/v1/vendor-credentials/`):

| Method | Path | Purpose |
|--------|------|---------|
| POST   | `/`            | save / update credential (email+password) |
| GET    | `/`            | list (masked) |
| DELETE | `/{id}`        | remove |

## Quickstart (UI flow)

1. **Save credential** (one-time): `POST /vendor-credentials` with
   `{vendor: "lc_bio", email, password, label: "default"}`.
2. **Open session**: `POST /data-downloads/sessions` with either
   `{vendor, email, password}` or `{vendor, credential_id}`. Returns
   `{active: true, pid: 12345}`.
3. **Find delivery on vendor's web UI**, click the linux-download button,
   copy the `obs_path` (the part inside the quotes after `download`).
4. **Start download**: `POST /data-downloads/jobs` with
   `{vendor: "lc_bio", source_path: "<obs_path>", dest_path: "/data/.../rawdata", auto_register: true, project_name?: "..."}`.
   Returns a `DownloadJob` with id + log_path.
5. **Watch progress**: subscribe via WebSocket to
   `realtime:task:<job_id>` (frontend uses the existing `/ws` endpoint
   with `subscribe_to_task` on the job's id), or poll
   `GET /data-downloads/jobs/{id}`.
6. **On completion**: if `auto_register=true`, the worker untars the
   delivery into `dest_path/extracted/` and creates an NGSmodule
   `Project` linked back to the download. The project's `config`
   carries `fastq_inventory` for the user to review and bind to samples.

## 联川 (lcbio) specifics

The Java daemon (`obs_service-...-SNAPSHOT.jar`) is **process-level
singleton** — only one session per backend host. Concurrent downloads
inside one session are sequential (the daemon reads `linuxDown.txt`
which is overwritten by each new request).

Configure path with `LC_BIO_BIN_DIR` (default
`/home/zhanghao/programs/lcbio/linux_download/bin`).

Log format parsed for progress (see `services/data_provider/lc_bio.py`):

```
已开始下载，请耐心等待......
<obs_path> ：下载中......       <ts>
<obs_path> ：已下载：50.00 MB  0.39%      <ts>     ← progress line
<obs_path> ：已下载：1.03 GB  8.15%       <ts>
<obs_path> ： 下载已完成！      <ts>                ← terminal marker
```

## Adding a new vendor

1. Subclass `DataProviderBase` in `app/services/data_provider/<name>.py`
   and implement `session_status / login / logout / start_download /
   parse_progress`.
2. Register in `app/services/data_provider/factory.py`.
3. The rest of the stack (DB, API, worker, WebSocket) requires no
   changes — `DataDownloadService` is vendor-agnostic.

A `NovogeneProvider` stub lives in `novogene.py`; the actual delivery
channel (sftp / OBS / vendor CLI) needs to be confirmed before that
adapter can be filled in.

## Security

- Vendor passwords saved via `POST /vendor-credentials` are
  **fernet-encrypted at rest** with a key derived from
  `settings.SECRET_KEY` (SHA-256 → base64). Rotating SECRET_KEY
  invalidates stored ciphertexts (intentional — forces re-entry).
- Plaintext passwords are **never returned by the API**; only masked
  email previews (`us***@example.com`).
- `linuxLogin.txt` (lcbio's plaintext credential file) is owned by the
  user that runs the backend, NOT by the API caller. The lcbio Java
  service truncates it to 0 bytes after reading, so it's only briefly
  on disk; restrict the backend host's filesystem accordingly.

## Known limitations

- lcbio's `download` CLI has no list-deliveries endpoint, so the
  frontend can't enumerate available files — the user must paste the
  `obs_path` copied from the vendor's web UI.
- `cancel` for lcbio kills the entire Java daemon (and any other
  in-flight downloads sharing the session). For finer control we'd
  need to call OBS directly, bypassing lcbio.
- Stale job timeout is hard-coded at 10 min without progress; tune
  `_STALE_TIMEOUT` in `workers/download_tasks.py` if you regularly
  download from slow OBS regions.

## Deployment notes

### lcbio container mount (required for 联川 downloads)

The lcbio bash CLI (`run` / `download` / `stop`) ships with an embedded
JRE and persists state under `linux_download/lib/.obs_service/`. The
backend container needs **read-write** access to this tree because the
Java daemon writes `linuxLogin.txt` / `linuxDown.txt` / `log.txt` and
expects to spawn the JVM from the bundled `jre1.8.0_431/`.

`docker-compose.yml` and `docker-compose.prod.yml` mount
`${LCBIO_HOST_DIR:-/opt/lcbio}:/opt/lcbio` and set
`LC_BIO_BIN_DIR=/opt/lcbio/linux_download/bin` for both the backend API
and the celery worker. Set `LCBIO_HOST_DIR` in your `.env` to wherever
lcbio is installed on the host (the official installer is
`https://www.omicstudio.cn/attached/file/newObs/lcbio_download.run`).

If you don't use 联川 deliveries, the default `/opt/lcbio` path is fine
even when nothing is there — the data-downloads endpoints will simply
return errors when invoked, but unrelated features keep working.

### Single-replica limitation

The lcbio Java daemon is **process-level singleton** (the `run` script
refuses to start a second instance via `ps grep`). Concurrent downloads
inside one session are sequential because `linuxDown.txt` is overwritten
on each request. **Don't scale the backend container beyond 1 replica**
until the daemon is moved to a dedicated sidecar service that all
backend replicas talk to.

To scale: extract the daemon into its own `lcbio-daemon` compose
service, share `/opt/lcbio/linux_download/lib/.obs_service/` between
that container and the backend via a named volume, and have backend
replicas write to `linuxDown.txt` with file locking.

### alembic migration chain

The repo has no baseline-create-tables migration; the existing chain
(`add_notifications_001` → `admin_enhanced_001` → `data_downloads_001`)
contains only incremental ALTERs. First-time deployments must therefore:

```bash
# 1. let SQLAlchemy create the canonical schema from models
python -c "from app.core.database import init_db; init_db()"

# 2. tell alembic the DB is already at HEAD (skips the chain)
alembic stamp head
```

Subsequent schema changes can use `alembic upgrade head` normally.

`alembic/env.py` reads `DATABASE_URL` from the environment when set,
so containerised deployments don't need to template `alembic.ini`.
