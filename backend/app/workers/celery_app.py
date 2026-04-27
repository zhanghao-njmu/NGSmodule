"""
Celery application configuration
"""
from celery import Celery
from celery.schedules import crontab
from app.core.config import settings

# Create Celery app
celery_app = Celery(
    "ngsmodule",
    broker=settings.CELERY_BROKER_URL,
    backend=settings.CELERY_RESULT_BACKEND,
)

# Configuration
celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    task_time_limit=86400,  # 24 hours
    task_soft_time_limit=82800,  # 23 hours
)

# Periodic tasks (Celery Beat)
celery_app.conf.beat_schedule = {
    "system-health-check": {
        "task": "admin.check_system_health",
        "schedule": 60.0,  # Every minute
    },
    "sync-celery-jobs": {
        "task": "admin.sync_celery_jobs",
        "schedule": 300.0,  # Every 5 minutes
    },
    "cleanup-expired-backups": {
        "task": "admin.cleanup_expired_backups",
        "schedule": crontab(hour=3, minute=0),  # Daily at 3 AM
    },
    "cleanup-old-alerts": {
        "task": "admin.cleanup_old_alerts",
        "schedule": crontab(hour=4, minute=0, day_of_week=0),  # Weekly on Sunday at 4 AM
    },
}

# Auto-discover tasks
celery_app.autodiscover_tasks(["app.workers"])
