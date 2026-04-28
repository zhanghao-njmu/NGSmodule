"""
Data download service - DB-side bookkeeping for vendor downloads.

The service owns the `DownloadJob` lifecycle. Vendor-specific work
(spawning the lcbio daemon, parsing logs) is delegated to a
`DataProviderBase` from `app.services.data_provider`.
"""
from pathlib import Path
from typing import List, Optional
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy.orm import Session

from app.core.datetime_utils import utc_now_naive
from app.models.data_download import DownloadJob
from app.services.data_provider import (
    get_provider,
    SessionExistsError,
    NoSessionError,
    LoginFailedError,
    DownloadFailedError,
)
from app.services.data_provider.base import SessionInfo


class DataDownloadService:
    """Owns the DownloadJob table and forwards vendor-side work to providers."""

    def __init__(self, db: Session):
        self.db = db

    # ---------------- session ----------------

    def session_status(self, vendor: str) -> SessionInfo:
        return get_provider(vendor).session_status()

    def login(self, vendor: str, email: str, password: str) -> SessionInfo:
        try:
            return get_provider(vendor).login(email, password)
        except SessionExistsError as exc:
            raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail=str(exc))
        except LoginFailedError as exc:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=str(exc))

    def logout(self, vendor: str) -> None:
        get_provider(vendor).logout()

    # ---------------- jobs ----------------

    def create_job(
        self,
        user_id: UUID,
        vendor: str,
        source_path: str,
        dest_path: str,
        auto_register: bool = False,
        project_name: Optional[str] = None,
    ) -> DownloadJob:
        """Persist a job and start the download. Returns the row.

        Worker tasks (Phase 2) will pick up persistence + tail parsing;
        for now we trigger the vendor synchronously and store the log path
        so progress polling works immediately.
        """
        provider = get_provider(vendor)
        try:
            log_path = provider.start_download(source_path, Path(dest_path))
        except NoSessionError as exc:
            raise HTTPException(status_code=status.HTTP_409_CONFLICT, detail=str(exc))
        except DownloadFailedError as exc:
            raise HTTPException(status_code=status.HTTP_502_BAD_GATEWAY, detail=str(exc))

        job = DownloadJob(
            user_id=user_id,
            vendor=vendor,
            source_path=source_path,
            dest_path=dest_path,
            log_path=str(log_path),
            status="running",
            started_at=utc_now_naive(),
            auto_register="true" if auto_register else "false",
            project_name_hint=project_name,
        )
        self.db.add(job)
        self.db.commit()
        self.db.refresh(job)

        # Hand off to Celery worker for long-running progress tailing +
        # realtime event publishing. Best-effort: if the broker is down
        # the row is still persisted and the API still works (clients can
        # poll GET /jobs/:id which calls refresh_progress synchronously).
        try:
            from app.workers.download_tasks import watch_progress
            watch_progress.delay(str(job.id))
        except Exception:
            # Worker enqueue failure is logged but doesn't fail the request.
            import logging
            logging.getLogger(__name__).warning(
                "Failed to enqueue watch_progress task; job will rely on sync polling",
                exc_info=True,
            )

        return job

    def get_by_id(self, job_id: UUID, user_id: UUID) -> Optional[DownloadJob]:
        return (
            self.db.query(DownloadJob)
            .filter(DownloadJob.id == job_id, DownloadJob.user_id == user_id)
            .first()
        )

    def get_by_id_or_raise(self, job_id: UUID, user_id: UUID) -> DownloadJob:
        job = self.get_by_id(job_id, user_id)
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Download job {job_id} not found",
            )
        return job

    def list_jobs(self, user_id: UUID) -> List[DownloadJob]:
        return (
            self.db.query(DownloadJob)
            .filter(DownloadJob.user_id == user_id)
            .order_by(DownloadJob.created_at.desc())
            .all()
        )

    def refresh_progress(self, job: DownloadJob) -> DownloadJob:
        """Re-parse the vendor log and update job state. Idempotent."""
        if job.status in ("completed", "failed", "cancelled") or not job.log_path:
            return job
        provider = get_provider(job.vendor)
        snap = provider.parse_progress(Path(job.log_path))

        job.progress_pct = snap.percent
        if snap.bytes_downloaded is not None:
            job.bytes_downloaded = snap.bytes_downloaded
        if snap.file_size is not None:
            job.file_size = snap.file_size
        if snap.status == "completed" and job.status != "completed":
            job.status = "completed"
            job.finished_at = utc_now_naive()
        elif snap.status == "failed" and job.status != "failed":
            job.status = "failed"
            job.error_message = snap.error_message
            job.finished_at = utc_now_naive()
        else:
            job.status = snap.status

        self.db.commit()
        self.db.refresh(job)
        return job

    def cancel(self, job: DownloadJob) -> DownloadJob:
        """Cancel a running job. Implementation note: lcbio's `stop` kills the
        whole daemon — caller should warn that this affects all in-flight
        downloads under the same session."""
        if job.status in ("completed", "failed", "cancelled"):
            return job
        get_provider(job.vendor).logout()
        job.status = "cancelled"
        job.finished_at = utc_now_naive()
        self.db.commit()
        self.db.refresh(job)
        return job
