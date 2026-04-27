"""
Alert Service
Real implementation for system alert monitoring and persistence.
"""
import logging
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta
from sqlalchemy.orm import Session
from sqlalchemy import desc, and_

from app.core.config import settings
from app.core.datetime_utils import utc_now_naive
from app.models.alert import SystemAlert


logger = logging.getLogger(__name__)


class AlertService:
    """Service for system alert operations"""

    def __init__(self, db: Session):
        self.db = db

    def create_alert(
        self,
        type: str,
        severity: str,
        title: str,
        message: str,
        source: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
        deduplicate_within_minutes: int = 15,
    ) -> SystemAlert:
        """
        Create a new alert with deduplication

        If a similar unresolved alert was created within the deduplication window,
        update the existing alert's last_seen and occurrence_count instead.
        """
        if deduplicate_within_minutes > 0:
            since = utc_now_naive() - timedelta(minutes=deduplicate_within_minutes)
            existing = (
                self.db.query(SystemAlert)
                .filter(
                    and_(
                        SystemAlert.title == title,
                        SystemAlert.source == source,
                        SystemAlert.resolved == False,
                        SystemAlert.last_seen >= since,
                    )
                )
                .first()
            )

            if existing:
                existing.last_seen = utc_now_naive()
                try:
                    existing.occurrence_count = str(int(existing.occurrence_count) + 1)
                except (ValueError, TypeError):
                    existing.occurrence_count = "2"
                existing.message = message  # Update with latest message
                self.db.commit()
                self.db.refresh(existing)
                return existing

        alert = SystemAlert(
            type=type,
            severity=severity,
            title=title,
            message=message,
            source=source,
            alert_metadata=metadata or {},
        )
        self.db.add(alert)
        self.db.commit()
        self.db.refresh(alert)
        return alert

    def get_alerts(
        self,
        resolved: Optional[bool] = False,
        type: Optional[str] = None,
        severity: Optional[str] = None,
        source: Optional[str] = None,
        skip: int = 0,
        limit: int = 100,
    ) -> List[SystemAlert]:
        """Query alerts with filters"""
        query = self.db.query(SystemAlert)

        if resolved is not None:
            query = query.filter(SystemAlert.resolved == resolved)
        if type:
            query = query.filter(SystemAlert.type == type)
        if severity:
            query = query.filter(SystemAlert.severity == severity)
        if source:
            query = query.filter(SystemAlert.source == source)

        return (
            query.order_by(desc(SystemAlert.timestamp))
            .offset(skip)
            .limit(limit)
            .all()
        )

    def count_alerts(self, resolved: Optional[bool] = None) -> int:
        """Count alerts matching filter"""
        query = self.db.query(SystemAlert)
        if resolved is not None:
            query = query.filter(SystemAlert.resolved == resolved)
        return query.count()

    def get_alert(self, alert_id: str) -> Optional[SystemAlert]:
        """Get alert by ID"""
        return self.db.query(SystemAlert).filter(SystemAlert.id == alert_id).first()

    def resolve_alert(
        self,
        alert_id: str,
        admin_user_id: str,
        notes: Optional[str] = None,
    ) -> Optional[SystemAlert]:
        """Mark an alert as resolved"""
        alert = self.get_alert(alert_id)
        if not alert:
            return None

        alert.resolved = True
        alert.resolved_at = utc_now_naive()
        alert.resolved_by = admin_user_id
        alert.resolution_notes = notes

        self.db.commit()
        self.db.refresh(alert)
        return alert

    def check_system_health(self) -> List[SystemAlert]:
        """
        Check system health metrics and generate alerts as needed.

        This should be run periodically (e.g., every minute via Celery beat).
        """
        import psutil
        new_alerts = []

        # Check disk space
        try:
            disk = psutil.disk_usage('/')
            if disk.percent >= settings.ALERT_DISK_CRITICAL_PERCENT:
                alert = self.create_alert(
                    type="error",
                    severity="critical",
                    title="Disk Space Critical",
                    message=f"Disk usage is at {disk.percent:.1f}%. Immediate action required.",
                    source="system_monitor",
                    metadata={
                        "disk_percent": disk.percent,
                        "free_gb": round(disk.free / (1024**3), 2),
                        "threshold": settings.ALERT_DISK_CRITICAL_PERCENT,
                    },
                )
                new_alerts.append(alert)
            elif disk.percent >= settings.ALERT_DISK_WARNING_PERCENT:
                alert = self.create_alert(
                    type="warning",
                    severity="high",
                    title="Disk Space Warning",
                    message=f"Disk usage is at {disk.percent:.1f}%. Consider cleanup.",
                    source="system_monitor",
                    metadata={
                        "disk_percent": disk.percent,
                        "free_gb": round(disk.free / (1024**3), 2),
                        "threshold": settings.ALERT_DISK_WARNING_PERCENT,
                    },
                )
                new_alerts.append(alert)
        except Exception as e:
            logger.error(f"Disk check failed: {e}")

        # Check memory
        try:
            memory = psutil.virtual_memory()
            if memory.percent >= settings.ALERT_MEMORY_CRITICAL_PERCENT:
                alert = self.create_alert(
                    type="error",
                    severity="critical",
                    title="Memory Critical",
                    message=f"Memory usage is at {memory.percent:.1f}%.",
                    source="system_monitor",
                    metadata={
                        "memory_percent": memory.percent,
                        "available_gb": round(memory.available / (1024**3), 2),
                        "threshold": settings.ALERT_MEMORY_CRITICAL_PERCENT,
                    },
                )
                new_alerts.append(alert)
            elif memory.percent >= settings.ALERT_MEMORY_WARNING_PERCENT:
                alert = self.create_alert(
                    type="warning",
                    severity="high",
                    title="High Memory Usage",
                    message=f"Memory usage is at {memory.percent:.1f}%.",
                    source="system_monitor",
                    metadata={
                        "memory_percent": memory.percent,
                        "available_gb": round(memory.available / (1024**3), 2),
                        "threshold": settings.ALERT_MEMORY_WARNING_PERCENT,
                    },
                )
                new_alerts.append(alert)
        except Exception as e:
            logger.error(f"Memory check failed: {e}")

        return new_alerts

    def cleanup_old_alerts(self, days_to_keep: int = 90) -> int:
        """Delete resolved alerts older than the specified number of days"""
        cutoff = utc_now_naive() - timedelta(days=days_to_keep)
        count = (
            self.db.query(SystemAlert)
            .filter(
                and_(
                    SystemAlert.resolved == True,
                    SystemAlert.resolved_at < cutoff,
                )
            )
            .delete(synchronize_session=False)
        )
        self.db.commit()
        return count
