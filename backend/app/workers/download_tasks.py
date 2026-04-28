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
    """Hook for Phase 5 (auto-extract + register as NGSmodule project).

    In Phase 2 this is a noop so the worker doesn't yet depend on
    project_service / sample_service. Phase 5 will replace it with
    untar + Project/Sample creation; signature is kept stable.
    """
    _ = (db, job)  # silence unused-warning until Phase 5
