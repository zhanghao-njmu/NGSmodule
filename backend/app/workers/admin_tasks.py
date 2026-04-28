"""
Admin Celery tasks for backups, monitoring, and cleanup
"""

import logging

from celery import shared_task
from sqlalchemy.orm import Session

from app.core.database import SessionLocal

logger = logging.getLogger(__name__)


def _get_db() -> Session:
    """Get a database session for tasks"""
    return SessionLocal()


@shared_task(bind=True, name="admin.create_backup")
def create_backup_task(
    self,
    backup_id: str,
    backup_type: str,
    admin_user_id: str,
    description: str = None,
    compress: bool = True,
):
    """
    Async backup creation task

    Updates the existing SystemBackup record (created with status='pending')
    by performing the backup and updating its final state.
    """
    from app.core.datetime_utils import utc_now_naive
    from app.models.backup import SystemBackup
    from app.services.backup_service import BackupService
    from app.services.job_service import JobService

    db = _get_db()
    try:
        # Update backup record to in_progress
        backup = db.query(SystemBackup).filter(SystemBackup.id == backup_id).first()
        if not backup:
            logger.error(f"Backup {backup_id} not found")
            return {"success": False, "error": "Backup record not found"}

        backup.status = "in_progress"
        backup.started_at = utc_now_naive()
        db.commit()

        job_service = JobService(db)

        # Perform actual backup
        backup_service = BackupService(db)
        result = backup_service.create_backup(
            backup_type=backup_type,
            admin_user_id=admin_user_id,
            description=description,
            compress=compress,
        )

        # Sync job status linked via celery_task_id
        celery_job = job_service.get_job_by_celery_id(self.request.id)
        if celery_job:
            if result.status == "completed":
                job_service.update_job_status(
                    str(celery_job.id),
                    status="completed",
                    progress=1.0,
                    result={"backup_id": str(result.id), "size": result.size},
                )
            else:
                job_service.update_job_status(
                    str(celery_job.id),
                    status="failed",
                    error=result.error_message or "Backup failed",
                )

        return {
            "success": result.status == "completed",
            "backup_id": str(result.id),
            "status": result.status,
            "size": result.size,
        }

    except Exception as e:
        logger.error(f"Backup task failed: {e}", exc_info=True)
        return {"success": False, "error": str(e)}
    finally:
        db.close()


@shared_task(name="admin.check_system_health")
def check_system_health_task():
    """
    Periodic system health check task.

    Runs every minute to check disk, memory, and other metrics.
    Generates alerts when thresholds are exceeded.
    """
    from app.services.alert_service import AlertService

    db = _get_db()
    try:
        service = AlertService(db)
        new_alerts = service.check_system_health()
        return {
            "alerts_created": len(new_alerts),
            "alert_ids": [str(a.id) for a in new_alerts],
        }
    except Exception as e:
        logger.error(f"Health check task failed: {e}", exc_info=True)
        return {"error": str(e)}
    finally:
        db.close()


@shared_task(name="admin.cleanup_expired_backups")
def cleanup_expired_backups_task():
    """
    Periodic task to clean up expired backups.

    Runs daily to delete backup files past their retention period.
    """
    from app.services.backup_service import BackupService

    db = _get_db()
    try:
        service = BackupService(db)
        count = service.cleanup_expired_backups()
        logger.info(f"Cleaned up {count} expired backups")
        return {"deleted_count": count}
    except Exception as e:
        logger.error(f"Backup cleanup task failed: {e}", exc_info=True)
        return {"error": str(e)}
    finally:
        db.close()


@shared_task(name="admin.cleanup_old_alerts")
def cleanup_old_alerts_task(days_to_keep: int = 90):
    """
    Periodic task to clean up old resolved alerts.

    Runs weekly to delete resolved alerts older than retention period.
    """
    from app.services.alert_service import AlertService

    db = _get_db()
    try:
        service = AlertService(db)
        count = service.cleanup_old_alerts(days_to_keep=days_to_keep)
        logger.info(f"Cleaned up {count} old alerts")
        return {"deleted_count": count}
    except Exception as e:
        logger.error(f"Alert cleanup task failed: {e}", exc_info=True)
        return {"error": str(e)}
    finally:
        db.close()


@shared_task(name="admin.sync_celery_jobs")
def sync_celery_jobs_task():
    """
    Periodic task to sync job statuses with Celery.

    Detects jobs that may have been lost or where state diverged.
    """
    from app.models.system_job import SystemJob
    from app.services.job_service import JobService

    db = _get_db()
    try:
        service = JobService(db)
        # Get all running/pending jobs
        active_jobs = (
            db.query(SystemJob)
            .filter(SystemJob.status.in_(["pending", "running"]))
            .filter(SystemJob.celery_task_id.isnot(None))
            .all()
        )

        synced_count = 0
        for job in active_jobs:
            service.sync_with_celery(job)
            synced_count += 1

        return {"synced_count": synced_count}
    except Exception as e:
        logger.error(f"Job sync task failed: {e}", exc_info=True)
        return {"error": str(e)}
    finally:
        db.close()
