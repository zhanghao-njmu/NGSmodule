"""
Job Service
Real implementation for system job tracking with Celery integration.
"""

import logging
from typing import Any, Dict, List, Optional

from sqlalchemy import desc
from sqlalchemy.orm import Session

from app.core.datetime_utils import utc_now_naive
from app.models.system_job import SystemJob

logger = logging.getLogger(__name__)


class JobService:
    """Service for system job tracking and management"""

    def __init__(self, db: Session):
        self.db = db

    def create_job(
        self,
        type: str,
        user_id: Optional[str] = None,
        username: Optional[str] = None,
        celery_task_id: Optional[str] = None,
        parameters: Optional[Dict[str, Any]] = None,
        message: Optional[str] = None,
        parent_job_id: Optional[str] = None,
    ) -> SystemJob:
        """Create a new system job record"""
        job = SystemJob(
            type=type,
            status="pending",
            user_id=user_id,
            username=username,
            celery_task_id=celery_task_id,
            parameters=parameters or {},
            message=message,
            parent_job_id=parent_job_id,
        )
        self.db.add(job)
        self.db.commit()
        self.db.refresh(job)
        return job

    def update_job_status(
        self,
        job_id: str,
        status: str,
        progress: Optional[float] = None,
        message: Optional[str] = None,
        result: Optional[Dict[str, Any]] = None,
        error: Optional[str] = None,
    ) -> Optional[SystemJob]:
        """Update job status and progress"""
        job = self.get_job(job_id)
        if not job:
            return None

        job.status = status

        if progress is not None:
            job.progress = max(0.0, min(1.0, progress))
        if message is not None:
            job.message = message
        if result is not None:
            job.result = result
        if error is not None:
            job.error = error

        # Update timestamps based on status
        now = utc_now_naive()
        if status == "running" and not job.started_at:
            job.started_at = now
        elif status in ("completed", "failed", "cancelled"):
            job.completed_at = now
            if status == "completed":
                job.progress = 1.0

        self.db.commit()
        self.db.refresh(job)
        return job

    def list_jobs(
        self,
        status: Optional[str] = None,
        type: Optional[str] = None,
        user_id: Optional[str] = None,
        skip: int = 0,
        limit: int = 50,
    ) -> List[SystemJob]:
        """List jobs with filters"""
        query = self.db.query(SystemJob)

        if status:
            query = query.filter(SystemJob.status == status)
        if type:
            query = query.filter(SystemJob.type == type)
        if user_id:
            query = query.filter(SystemJob.user_id == user_id)

        return query.order_by(desc(SystemJob.created_at)).offset(skip).limit(limit).all()

    def count_jobs(
        self,
        status: Optional[str] = None,
        type: Optional[str] = None,
    ) -> int:
        """Count jobs matching filters"""
        query = self.db.query(SystemJob)
        if status:
            query = query.filter(SystemJob.status == status)
        if type:
            query = query.filter(SystemJob.type == type)
        return query.count()

    def get_job(self, job_id: str) -> Optional[SystemJob]:
        """Get job by ID"""
        return self.db.query(SystemJob).filter(SystemJob.id == job_id).first()

    def get_job_by_celery_id(self, celery_task_id: str) -> Optional[SystemJob]:
        """Get job by Celery task ID"""
        return self.db.query(SystemJob).filter(SystemJob.celery_task_id == celery_task_id).first()

    def cancel_job(self, job_id: str) -> Optional[SystemJob]:
        """
        Cancel a running or pending job.

        Revokes the Celery task if applicable and updates the job status.
        """
        job = self.get_job(job_id)
        if not job:
            return None

        if job.status not in ("pending", "running"):
            logger.warning(f"Job {job_id} is in {job.status} state, cannot cancel")
            return job

        # Revoke Celery task if exists
        if job.celery_task_id:
            try:
                from app.workers.celery_app import celery_app

                celery_app.control.revoke(job.celery_task_id, terminate=True)
                logger.info(f"Revoked Celery task {job.celery_task_id}")
            except Exception as e:
                logger.error(f"Failed to revoke Celery task: {e}")

        job.status = "cancelled"
        job.completed_at = utc_now_naive()
        job.message = "Job cancelled by administrator"
        self.db.commit()
        self.db.refresh(job)
        return job

    def retry_job(self, job_id: str) -> Optional[SystemJob]:
        """
        Retry a failed job by creating a new job with the same parameters.
        """
        original = self.get_job(job_id)
        if not original:
            return None

        if original.status != "failed":
            logger.warning(f"Job {job_id} is in {original.status} state, only failed jobs can be retried")
            return None

        try:
            new_retry_count = str(int(original.retry_count) + 1)
        except (ValueError, TypeError):
            new_retry_count = "1"

        new_job = SystemJob(
            type=original.type,
            status="pending",
            user_id=original.user_id,
            username=original.username,
            parameters=original.parameters,
            message=f"Retry of job {original.id}",
            retry_count=new_retry_count,
            parent_job_id=original.id,
        )
        self.db.add(new_job)
        self.db.commit()
        self.db.refresh(new_job)
        return new_job

    def sync_with_celery(self, job: SystemJob) -> SystemJob:
        """
        Sync job status with Celery state.

        Useful for jobs that were lost during a crash or where the
        database state doesn't match the actual Celery state.
        """
        if not job.celery_task_id:
            return job

        try:
            from celery.result import AsyncResult

            from app.workers.celery_app import celery_app

            result = AsyncResult(job.celery_task_id, app=celery_app)
            celery_state = result.state

            # Map Celery states to our job states
            state_map = {
                "PENDING": "pending",
                "STARTED": "running",
                "SUCCESS": "completed",
                "FAILURE": "failed",
                "RETRY": "running",
                "REVOKED": "cancelled",
            }

            new_status = state_map.get(celery_state)
            if new_status and new_status != job.status:
                job.status = new_status

                if celery_state == "SUCCESS":
                    job.progress = 1.0
                    job.completed_at = utc_now_naive()
                    if result.result is not None:
                        try:
                            job.result = (
                                result.result if isinstance(result.result, dict) else {"output": str(result.result)}
                            )
                        except Exception:
                            pass
                elif celery_state == "FAILURE":
                    job.completed_at = utc_now_naive()
                    job.error = str(result.info) if result.info else "Task failed"
                elif celery_state == "STARTED" and not job.started_at:
                    job.started_at = utc_now_naive()

                self.db.commit()
                self.db.refresh(job)
        except Exception as e:
            logger.error(f"Failed to sync job {job.id} with Celery: {e}")

        return job

    def get_active_jobs_count(self) -> int:
        """Count jobs that are pending or running"""
        return self.db.query(SystemJob).filter(SystemJob.status.in_(["pending", "running"])).count()
