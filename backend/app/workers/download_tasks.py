"""
Celery tasks for vendor data downloads.

The actual byte-shuffling happens inside the vendor's daemon (lcbio's
Java service); our task just tails the vendor log file, updates the
DownloadJob row, and publishes realtime progress events.
"""
import logging
import time
from pathlib import Path
from typing import Optional

from app.workers.celery_app import celery_app
from app.core.database import SessionLocal
from app.core.datetime_utils import utc_now_naive
from app.models.data_download import DownloadJob
from app.services.data_provider import get_provider

logger = logging.getLogger(__name__)


# How often to re-parse the vendor log. lcbio writes progress every ~5s
# at 50MB chunks, so 2s polling gives a smooth UI without being wasteful.
_POLL_INTERVAL = 2.0
# Hard cap to prevent runaway tasks if the vendor log stops updating but
# the daemon also doesn't write 'completed' (e.g. machine power-cycle).
_STALE_TIMEOUT = 600.0  # 10 min without progress change ⇒ mark failed
# Total task time limit (Celery worker enforces ~24h via celery_app config).
_MAX_RUNTIME = 24 * 3600.0


def _publish(job_id: str, status_str: str, progress: float, message: str = "") -> None:
    """Emit realtime event so connected WebSocket clients see progress.

    Reuses the per-task channel (`realtime:task:<id>`) which the existing
    websocket subscriber already routes; no new channel needed.
    """
    try:
        from app.services.realtime import publish_task_event

        publish_task_event(
            job_id,
            "download_update",
            {
                "job_id": job_id,
                "status": status_str,
                "progress": progress,
                "message": message,
            },
        )
    except Exception as exc:
        logger.debug(f"realtime publish skipped: {exc}")


@celery_app.task(bind=True, name="data_downloads.watch_progress")
def watch_progress(self, job_id: str) -> dict:
    """Tail the vendor log for `job_id` until the download terminates.

    Idempotent — safe to retry. If the job row is already in a terminal
    state, we noop.
    """
    db = SessionLocal()
    try:
        job: Optional[DownloadJob] = (
            db.query(DownloadJob).filter(DownloadJob.id == job_id).first()
        )
        if job is None:
            logger.warning(f"watch_progress: job {job_id} not found")
            return {"job_id": job_id, "status": "not_found"}
        if job.status in ("completed", "failed", "cancelled"):
            return {"job_id": job_id, "status": job.status}
        if not job.log_path:
            job.status = "failed"
            job.error_message = "missing log_path"
            job.finished_at = utc_now_naive()
            db.commit()
            return {"job_id": job_id, "status": "failed"}

        provider = get_provider(job.vendor)
        log_path = Path(job.log_path)
        start_time = time.time()
        last_change_at = start_time

        while True:
            if time.time() - start_time > _MAX_RUNTIME:
                job.status = "failed"
                job.error_message = f"exceeded max runtime ({_MAX_RUNTIME}s)"
                job.finished_at = utc_now_naive()
                db.commit()
                _publish(job_id, "failed", job.progress_pct, job.error_message)
                return {"job_id": job_id, "status": "failed"}

            snap = provider.parse_progress(log_path)

            # Persist any state change.
            changed = False
            if snap.percent != job.progress_pct:
                job.progress_pct = snap.percent
                changed = True
            if snap.bytes_downloaded is not None and snap.bytes_downloaded != job.bytes_downloaded:
                job.bytes_downloaded = snap.bytes_downloaded
                changed = True
            if snap.file_size is not None and snap.file_size != job.file_size:
                job.file_size = snap.file_size
                changed = True

            if snap.status == "completed":
                job.status = "completed"
                job.progress_pct = 100.0
                job.finished_at = utc_now_naive()
                db.commit()
                _publish(job_id, "completed", 100.0, "")
                _maybe_post_download_hook(db, job)
                return {"job_id": job_id, "status": "completed"}

            if snap.status == "failed":
                job.status = "failed"
                job.error_message = snap.error_message or "vendor reported failure"
                job.finished_at = utc_now_naive()
                db.commit()
                _publish(job_id, "failed", job.progress_pct, job.error_message)
                return {"job_id": job_id, "status": "failed"}

            if changed:
                last_change_at = time.time()
                job.status = snap.status  # "running" mostly
                db.commit()
                _publish(job_id, snap.status, snap.percent, "")
            elif time.time() - last_change_at > _STALE_TIMEOUT:
                job.status = "failed"
                job.error_message = (
                    f"no progress for {int(_STALE_TIMEOUT)}s; vendor daemon may have died"
                )
                job.finished_at = utc_now_naive()
                db.commit()
                _publish(job_id, "failed", job.progress_pct, job.error_message)
                return {"job_id": job_id, "status": "failed"}

            time.sleep(_POLL_INTERVAL)

            # Periodically re-fetch the row in case another process cancelled.
            db.refresh(job)
            if job.status == "cancelled":
                _publish(job_id, "cancelled", job.progress_pct, "")
                return {"job_id": job_id, "status": "cancelled"}

    finally:
        db.close()


def _maybe_post_download_hook(db, job: DownloadJob) -> None:
    """If `auto_register` was set when the job was created, untar the
    delivery archive into `dest_path/extracted/` and create an NGSmodule
    Project pointing at it. Samples are NOT auto-created — naming is
    vendor-specific, so we leave that for the user to review in the UI.

    Best-effort: any failure is logged + recorded in error_message but
    doesn't change the job's terminal status (it's already "completed"
    by the caller before this hook runs).
    """
    if (job.auto_register or "").lower() != "true":
        return

    import tarfile
    from pathlib import Path
    from app.schemas.project import ProjectCreate
    from app.services.project_service import ProjectService

    dest = Path(job.dest_path)
    archive = _find_archive(dest, basename=job.source_path.rstrip("/").split("/")[-1])
    if archive is None:
        logger.info(f"job {job.id}: no tar archive found in {dest}; skipping registration")
        return

    extract_dir = dest / "extracted"
    extract_dir.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive, "r:*") as tf:
            # Refuse anything that escapes the target dir (path traversal guard)
            for m in tf.getmembers():
                if m.name.startswith("/") or ".." in Path(m.name).parts:
                    raise ValueError(f"unsafe archive entry: {m.name!r}")
            tf.extractall(extract_dir)
    except (tarfile.TarError, ValueError, OSError) as exc:
        msg = f"auto-extract failed: {exc}"
        logger.warning(f"job {job.id}: {msg}")
        job.error_message = (job.error_message or "") + f" | {msg}"
        db.commit()
        return

    # Pick a project name: user override → first segment of source_path → archive stem
    name = (
        job.project_name_hint
        or job.source_path.split("/", 1)[0]
        or archive.stem
    )[:100]

    try:
        svc = ProjectService(db)
        project = svc.create(
            user_id=job.user_id,
            project_data=ProjectCreate(
                name=name,
                description=f"Auto-created from {job.vendor} download {job.source_path!r}",
                project_type=None,  # user picks type later via UI
                config={
                    "source": job.vendor,
                    "source_path": job.source_path,
                    "extracted_path": str(extract_dir),
                    "fastq_inventory": _inventory_fastqs(extract_dir),
                    "auto_registered_from_download": str(job.id),
                },
            ),
        )
        job.project_id = project.id
        db.commit()
        logger.info(f"job {job.id}: registered project {project.id} ({project.name})")
    except Exception as exc:
        # ProjectService raises HTTPException for validation errors; we
        # surface them to the user via error_message but don't fail the
        # job (download itself succeeded).
        msg = f"project auto-register failed: {exc}"
        logger.warning(f"job {job.id}: {msg}")
        job.error_message = (job.error_message or "") + f" | {msg}"
        db.commit()


def _find_archive(dest: Path, basename: str) -> Optional[Path]:
    """Return the most-recently-modified tar archive matching the delivery."""
    # Prefer exact basename (e.g. Data.tar) if present, otherwise any tar in dest
    exact = dest / basename
    if exact.exists() and exact.is_file():
        return exact
    candidates = list(dest.glob("*.tar")) + list(dest.glob("*.tar.gz")) + list(dest.glob("*.tgz"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def _inventory_fastqs(root: Path, limit: int = 1000) -> list[str]:
    """Collect fastq filenames so the user can review samples in the UI.

    Capped at `limit` to keep the JSONB column small even on jumbo deliveries.
    """
    patterns = ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq")
    found: list[str] = []
    for pat in patterns:
        for p in root.rglob(pat):
            found.append(str(p.relative_to(root)))
            if len(found) >= limit:
                return found
    return found
