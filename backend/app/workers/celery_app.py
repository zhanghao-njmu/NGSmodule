"""
Celery application configuration
"""

from celery import Celery
from celery.schedules import crontab
from kombu import Queue

from app.core.config import settings

# Create Celery app
celery_app = Celery(
    "ngsmodule",
    broker=settings.CELERY_BROKER_URL,
    backend=settings.CELERY_RESULT_BACKEND,
)

# Multi-queue configuration:
#   default       — short, latency-sensitive tasks (notifications, audit)
#   admin         — admin operations (cleanup, sync, alert checks)
#   backup        — backup creation (long-running)
#   ngs_pipeline  — NGS pipeline execution (very long-running, heavy)
celery_app.conf.task_queues = (
    Queue("default", routing_key="default"),
    Queue("admin", routing_key="admin"),
    Queue("backup", routing_key="backup"),
    Queue("ngs_pipeline", routing_key="ngs_pipeline"),
)

celery_app.conf.task_default_queue = "default"
celery_app.conf.task_default_exchange = "default"
celery_app.conf.task_default_routing_key = "default"

celery_app.conf.task_routes = {
    # Admin housekeeping
    "admin.create_backup": {"queue": "backup"},
    "admin.cleanup_expired_backups": {"queue": "admin"},
    "admin.cleanup_old_alerts": {"queue": "admin"},
    "admin.check_system_health": {"queue": "admin"},
    "admin.sync_celery_jobs": {"queue": "admin"},
    # NGS pipeline tasks
    "app.workers.pipeline_tasks.*": {"queue": "ngs_pipeline"},
    # Vendor data downloads (long-running tail loops, IO-bound)
    "data_downloads.*": {"queue": "default"},
}

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
    # Worker hardening for multi-replica deployments
    worker_prefetch_multiplier=1,  # fetch one task at a time so heavy tasks
    # don't get hoarded by a single worker
    task_acks_late=True,  # ack only after success → safer with replicas
    task_reject_on_worker_lost=True,
    broker_connection_retry_on_startup=True,
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
