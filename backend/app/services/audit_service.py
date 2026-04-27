"""
Audit Log Service
Provides centralized audit logging for sensitive operations.
"""
import logging
from typing import Optional, List, Dict, Any
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import desc

from app.models.audit_log import AuditLog


logger = logging.getLogger(__name__)


class AuditService:
    """Service for audit log operations"""

    def __init__(self, db: Session):
        self.db = db

    def log_action(
        self,
        action: str,
        admin_user_id: str,
        admin_username: str,
        target_user_id: Optional[str] = None,
        target_username: Optional[str] = None,
        target_resource_type: Optional[str] = None,
        target_resource_id: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
        status: str = "success",
        ip_address: Optional[str] = None,
        user_agent: Optional[str] = None,
        request_id: Optional[str] = None,
    ) -> AuditLog:
        """Record an audit log entry"""
        try:
            entry = AuditLog(
                action=action,
                admin_user_id=admin_user_id,
                admin_username=admin_username,
                target_user_id=target_user_id,
                target_username=target_username,
                target_resource_type=target_resource_type,
                target_resource_id=target_resource_id,
                details=details or {},
                status=status,
                ip_address=ip_address,
                user_agent=user_agent,
                request_id=request_id,
            )
            self.db.add(entry)
            self.db.commit()
            self.db.refresh(entry)
            return entry
        except Exception as e:
            logger.error(f"Failed to log audit action {action}: {e}")
            self.db.rollback()
            raise

    def get_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        admin_user_id: Optional[str] = None,
        target_user_id: Optional[str] = None,
        action: Optional[str] = None,
        status: Optional[str] = None,
        skip: int = 0,
        limit: int = 50,
    ) -> List[AuditLog]:
        """Query audit logs with filters"""
        query = self.db.query(AuditLog)

        if start_date:
            query = query.filter(AuditLog.timestamp >= start_date)
        if end_date:
            query = query.filter(AuditLog.timestamp <= end_date)
        if admin_user_id:
            query = query.filter(AuditLog.admin_user_id == admin_user_id)
        if target_user_id:
            query = query.filter(AuditLog.target_user_id == target_user_id)
        if action:
            query = query.filter(AuditLog.action == action)
        if status:
            query = query.filter(AuditLog.status == status)

        return (
            query.order_by(desc(AuditLog.timestamp))
            .offset(skip)
            .limit(limit)
            .all()
        )

    def count_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        admin_user_id: Optional[str] = None,
        action: Optional[str] = None,
    ) -> int:
        """Count audit logs matching filters"""
        query = self.db.query(AuditLog)

        if start_date:
            query = query.filter(AuditLog.timestamp >= start_date)
        if end_date:
            query = query.filter(AuditLog.timestamp <= end_date)
        if admin_user_id:
            query = query.filter(AuditLog.admin_user_id == admin_user_id)
        if action:
            query = query.filter(AuditLog.action == action)

        return query.count()

    def export_logs_to_csv(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        admin_user_id: Optional[str] = None,
        action: Optional[str] = None,
    ) -> str:
        """Export audit logs to CSV format string"""
        import csv
        import io

        logs = self.get_logs(
            start_date=start_date,
            end_date=end_date,
            admin_user_id=admin_user_id,
            action=action,
            skip=0,
            limit=100000,  # Large limit for export
        )

        output = io.StringIO()
        writer = csv.writer(output)

        # Header
        writer.writerow([
            "ID", "Timestamp", "Action", "Admin User ID", "Admin Username",
            "Target User ID", "Target Username", "Resource Type", "Resource ID",
            "Status", "IP Address", "Details"
        ])

        # Rows
        for log in logs:
            writer.writerow([
                str(log.id),
                log.timestamp.isoformat() if log.timestamp else "",
                log.action,
                str(log.admin_user_id),
                log.admin_username,
                str(log.target_user_id) if log.target_user_id else "",
                log.target_username or "",
                log.target_resource_type or "",
                log.target_resource_id or "",
                log.status,
                log.ip_address or "",
                str(log.details or {}),
            ])

        return output.getvalue()
