"""
Admin service for user management, system configuration, and logs
"""

import json
import logging
import math
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional

from sqlalchemy import and_, desc, func, or_
from sqlalchemy.orm import Session

from app.core.security import get_password_hash
from app.models.file import File
from app.models.notification import Notification
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.models.user import User
from app.schemas.admin import *

logger = logging.getLogger(__name__)


class AdminService:
    """Service for admin operations"""

    def __init__(self, db: Session):
        self.db = db

    # ========================================================================
    # User Management
    # ========================================================================

    def get_users(
        self,
        skip: int = 0,
        limit: int = 50,
        role: Optional[UserRole] = None,
        is_active: Optional[bool] = None,
        search: Optional[str] = None,
        sort_by: str = "created_at",
        sort_order: str = "desc",
    ) -> UserListResponse:
        """Get paginated list of all users"""

        query = self.db.query(User)

        # Apply filters
        if role:
            query = query.filter(User.role == role.value)

        if is_active is not None:
            query = query.filter(User.is_active == is_active)

        if search:
            search_term = f"%{search}%"
            query = query.filter(
                or_(
                    User.username.ilike(search_term),
                    User.email.ilike(search_term),
                    User.full_name.ilike(search_term),
                    User.organization.ilike(search_term),
                )
            )

        # Get total count
        total = query.count()

        # Apply sorting
        sort_column = getattr(User, sort_by, User.created_at)
        if sort_order == "desc":
            query = query.order_by(desc(sort_column))
        else:
            query = query.order_by(sort_column)

        # Apply pagination
        users = query.offset(skip).limit(limit).all()

        # Convert to response model
        user_list = []
        for user in users:
            user_list.append(
                AdminUserList(
                    id=str(user.id),
                    username=user.username,
                    email=user.email,
                    full_name=user.full_name,
                    role=user.role,
                    organization=user.organization,
                    is_active=user.is_active,
                    storage_used=user.storage_used,
                    storage_quota=user.storage_quota,
                    storage_percent=user.storage_percent_used,
                    last_login=user.last_login,
                    created_at=user.created_at,
                    updated_at=user.updated_at,
                )
            )

        total_pages = math.ceil(total / limit) if limit > 0 else 0
        page = (skip // limit) + 1 if limit > 0 else 1

        return UserListResponse(users=user_list, total=total, page=page, page_size=limit, total_pages=total_pages)

    def get_user_detail(self, user_id: str) -> Optional[AdminUserDetail]:
        """Get detailed information about a user"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        # Count user's resources
        total_projects = self.db.query(func.count(Project.id)).filter(Project.owner_id == user_id).scalar() or 0

        total_samples = (
            self.db.query(func.count(Sample.id)).join(Project).filter(Project.owner_id == user_id).scalar() or 0
        )

        total_tasks = (
            self.db.query(func.count(PipelineTask.id))
            .join(Sample)
            .join(Project)
            .filter(Project.owner_id == user_id)
            .scalar()
            or 0
        )

        total_files = (
            self.db.query(func.count(File.id)).join(Sample).join(Project).filter(Project.owner_id == user_id).scalar()
            or 0
        )

        # Get last activity
        last_activity = (
            self.db.query(func.max(PipelineTask.updated_at))
            .join(Sample)
            .join(Project)
            .filter(Project.owner_id == user_id)
            .scalar()
        )

        return AdminUserDetail(
            id=str(user.id),
            username=user.username,
            email=user.email,
            full_name=user.full_name,
            role=user.role,
            organization=user.organization,
            is_active=user.is_active,
            storage_used=user.storage_used,
            storage_quota=user.storage_quota,
            storage_percent=user.storage_percent_used,
            last_login=user.last_login,
            created_at=user.created_at,
            updated_at=user.updated_at,
            total_projects=total_projects,
            total_samples=total_samples,
            total_tasks=total_tasks,
            total_files=total_files,
            last_activity=last_activity,
            login_count=0,  # Would need login tracking table
        )

    def update_user(self, user_id: str, update_data: UserUpdateRequest) -> Optional[AdminUserDetail]:
        """Update user information"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        # Update fields
        update_dict = update_data.dict(exclude_unset=True)
        for field, value in update_dict.items():
            if hasattr(user, field) and value is not None:
                setattr(user, field, value)

        user.updated_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(user)

        return self.get_user_detail(str(user.id))

    def change_user_role(self, user_id: str, new_role: UserRole) -> Optional[AdminUserDetail]:
        """Change user role"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        user.role = new_role.value
        user.updated_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(user)

        logger.info(f"User {user.username} role changed to {new_role.value}")

        return self.get_user_detail(str(user.id))

    def activate_deactivate_user(
        self, user_id: str, is_active: bool, reason: Optional[str] = None
    ) -> Optional[AdminUserDetail]:
        """Activate or deactivate a user"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        user.is_active = is_active
        user.updated_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(user)

        action = "activated" if is_active else "deactivated"
        logger.info(f"User {user.username} {action}. Reason: {reason}")

        return self.get_user_detail(str(user.id))

    def reset_user_password(self, user_id: str, new_password: str, notify_user: bool = True) -> bool:
        """Reset user password"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False

        user.password_hash = get_password_hash(new_password)
        user.updated_at = datetime.utcnow()

        self.db.commit()

        logger.info(f"Password reset for user {user.username}")

        # Note: Email notification feature is not yet implemented
        # When implemented, this should send notification if notify_user is True

        return True

    def delete_user(self, user_id: str, transfer_data_to: Optional[str] = None) -> bool:
        """Delete a user (with optional data transfer)"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False

        # Transfer data if requested
        if transfer_data_to:
            target_user = self.db.query(User).filter(User.id == transfer_data_to).first()
            if target_user:
                # Transfer projects
                self.db.query(Project).filter(Project.owner_id == user_id).update({"owner_id": transfer_data_to})
                logger.info(f"Transferred user {user.username} data to {target_user.username}")

        # Delete user (cascade will handle related records)
        self.db.delete(user)
        self.db.commit()

        logger.info(f"User {user.username} deleted")

        return True

    # ========================================================================
    # System Configuration
    # ========================================================================

    def get_system_config(self) -> SystemConfig:
        """Get complete system configuration"""

        # In production, this would read from a config table or file
        # For now, return default values

        return SystemConfig(
            general={
                "app_name": "NGSmodule",
                "app_version": "1.0.0",
                "maintenance_mode": False,
                "allow_registration": True,
                "default_user_quota": 107374182400,  # 100GB
                "max_upload_size": 524288000,  # 500MB
                "session_timeout": 3600,  # 1 hour
            },
            security={
                "jwt_expiry": 1800,  # 30 minutes
                "refresh_token_expiry": 604800,  # 7 days
                "password_min_length": 8,
                "require_email_verification": False,
                "enable_2fa": False,
                "max_login_attempts": 5,
                "lockout_duration": 900,  # 15 minutes
            },
            storage={
                "storage_backend": "minio",
                "storage_path": "/data/storage",
                "enable_compression": True,
                "enable_deduplication": False,
                "cleanup_threshold": 90,  # percentage
                "retention_policy_days": 365,
            },
            email={
                "smtp_host": "smtp.gmail.com",
                "smtp_port": 587,
                "smtp_use_tls": True,
                "smtp_from": "noreply@ngsmodule.com",
                "enable_email_notifications": True,
            },
            notification={
                "enable_notifications": True,
                "default_email_enabled": True,
                "default_app_enabled": True,
                "default_push_enabled": False,
                "notification_retention_days": 90,
            },
            pipeline={
                "max_concurrent_jobs": 4,
                "job_timeout": 86400,  # 24 hours
                "enable_auto_retry": True,
                "max_retries": 3,
                "default_pipeline": "rna_seq",
            },
            performance={
                "api_rate_limit": 100,  # requests per minute
                "db_pool_size": 20,
                "db_max_overflow": 10,
                "cache_enabled": True,
                "cache_ttl": 300,  # 5 minutes
            },
            last_updated=datetime.utcnow(),
        )

    def update_system_config(
        self, category: SystemConfigCategory, updates: Dict[str, Any], admin_id: str
    ) -> SystemConfig:
        """Update system configuration"""

        # In production, this would update a config table or file
        logger.info(f"System config {category.value} updated by admin {admin_id}: {updates}")

        # Note: Config persistence using database or file storage is not yet implemented
        # Currently returns in-memory defaults; updates are logged but not persisted

        return self.get_system_config()

    def reset_system_config(self, categories: Optional[List[SystemConfigCategory]] = None) -> SystemConfig:
        """Reset system configuration to defaults"""

        if categories:
            logger.info(f"Resetting config categories: {[c.value for c in categories]}")
        else:
            logger.info("Resetting all config to defaults")

        # Note: Config reset functionality returns defaults; not yet connected to persistence layer

        return self.get_system_config()

    # ========================================================================
    # System Logs
    # ========================================================================

    def get_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        levels: Optional[List[LogLevel]] = None,
        sources: Optional[List[LogSource]] = None,
        search: Optional[str] = None,
        limit: int = 100,
        offset: int = 0,
    ) -> LogResponse:
        """Query system logs"""

        # In production, this would query from a log aggregation system
        # For now, read from log files

        logs = []

        # Mock log entries
        if not start_date:
            start_date = datetime.utcnow() - timedelta(days=1)
        if not end_date:
            end_date = datetime.utcnow()

        # Generate some sample logs
        log_levels = levels if levels else [LogLevel.INFO, LogLevel.WARNING, LogLevel.ERROR]
        log_sources = sources if sources else [LogSource.API, LogSource.DATABASE, LogSource.PIPELINE]

        sample_messages = [
            "API request completed successfully",
            "Database connection established",
            "Pipeline task started",
            "User authentication successful",
            "File upload completed",
            "Warning: High memory usage detected",
            "Error: Failed to connect to external service",
        ]

        # Create 10 sample log entries
        for i in range(min(10, limit)):
            logs.append(
                LogEntry(
                    timestamp=datetime.utcnow() - timedelta(minutes=i * 10),
                    level=log_levels[i % len(log_levels)],
                    source=log_sources[i % len(log_sources)],
                    message=sample_messages[i % len(sample_messages)],
                    details={"request_id": f"req_{i}", "duration_ms": 150 + i * 10},
                    ip_address="192.168.1.1",
                )
            )

        # Apply search filter
        if search:
            logs = [log for log in logs if search.lower() in log.message.lower()]

        total = len(logs)
        has_more = total > (offset + limit)

        # Apply pagination
        logs = logs[offset : offset + limit]

        return LogResponse(logs=logs, total=total, has_more=has_more)

    def download_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        levels: Optional[List[LogLevel]] = None,
        sources: Optional[List[LogSource]] = None,
        format: str = "json",
    ) -> str:
        """Generate log file for download"""

        logs_response = self.get_logs(
            start_date=start_date, end_date=end_date, levels=levels, sources=sources, limit=10000
        )

        # In production, this would create a file and return a download URL
        if format == "json":
            log_data = json.dumps([log.dict() for log in logs_response.logs], default=str, indent=2)
        elif format == "csv":
            # Simple CSV format
            log_data = "timestamp,level,source,message\n"
            for log in logs_response.logs:
                log_data += f"{log.timestamp},{log.level.value},{log.source.value},{log.message}\n"
        else:  # txt
            log_data = "\n".join(
                [
                    f"[{log.timestamp}] {log.level.value.upper()} [{log.source.value}] {log.message}"
                    for log in logs_response.logs
                ]
            )

        # Save to temp file and return path
        temp_file = f"/tmp/logs_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.{format}"
        with open(temp_file, "w") as f:
            f.write(log_data)

        return temp_file

    # ========================================================================
    # System Health & Maintenance
    # ========================================================================

    def get_system_health(self) -> SystemHealth:
        """Get system health status"""

        services = []

        # Check database
        try:
            self.db.execute("SELECT 1")
            services.append(
                ServiceHealth(
                    name="PostgreSQL",
                    status=ServiceStatus.HEALTHY,
                    response_time=5.2,
                    last_check=datetime.utcnow(),
                    message="Database is responsive",
                )
            )
        except Exception as e:
            services.append(
                ServiceHealth(
                    name="PostgreSQL",
                    status=ServiceStatus.DOWN,
                    last_check=datetime.utcnow(),
                    message=f"Database error: {str(e)}",
                )
            )

        # Check Redis (mock)
        services.append(
            ServiceHealth(
                name="Redis",
                status=ServiceStatus.HEALTHY,
                response_time=1.5,
                last_check=datetime.utcnow(),
                message="Cache is operational",
            )
        )

        # Check MinIO (mock)
        services.append(
            ServiceHealth(
                name="MinIO",
                status=ServiceStatus.HEALTHY,
                response_time=8.3,
                last_check=datetime.utcnow(),
                message="Object storage is available",
            )
        )

        # Check Celery (mock)
        services.append(
            ServiceHealth(
                name="Celery",
                status=ServiceStatus.HEALTHY,
                response_time=3.1,
                last_check=datetime.utcnow(),
                message="Task queue is processing",
            )
        )

        # Determine overall status
        if all(s.status == ServiceStatus.HEALTHY for s in services):
            overall_status = ServiceStatus.HEALTHY
        elif any(s.status == ServiceStatus.DOWN for s in services):
            overall_status = ServiceStatus.DEGRADED
        else:
            overall_status = ServiceStatus.HEALTHY

        return SystemHealth(
            status=overall_status,
            services=services,
            timestamp=datetime.utcnow(),
            uptime=86400.0,  # Mock 1 day uptime
            version="1.0.0",
        )

    def cleanup_system(self, options: CleanupOptions) -> CleanupResponse:
        """Perform system cleanup"""

        results = []
        total_items = 0
        total_space = 0
        start_time = datetime.utcnow()

        cutoff_date = datetime.utcnow() - timedelta(days=options.days_to_keep)

        # Clean old logs
        if options.old_logs:
            # In production, this would actually delete log files
            items_deleted = 0
            space_freed = 0

            if not options.dry_run:
                # Delete logs older than cutoff_date
                logger.info(f"Deleting logs older than {cutoff_date}")
                # Implementation would go here

            results.append(
                CleanupResult(
                    operation="old_logs", items_deleted=items_deleted, space_freed=space_freed, duration=0.5, errors=[]
                )
            )
            total_items += items_deleted
            total_space += space_freed

        # Clean temp files
        if options.temp_files:
            items_deleted = 0
            space_freed = 0

            if not options.dry_run:
                # Delete temp files
                temp_path = Path("/tmp")
                if temp_path.exists():
                    for file in temp_path.glob("ngsmodule_temp_*"):
                        if file.stat().st_mtime < cutoff_date.timestamp():
                            size = file.stat().st_size
                            file.unlink()
                            items_deleted += 1
                            space_freed += size

            results.append(
                CleanupResult(
                    operation="temp_files",
                    items_deleted=items_deleted,
                    space_freed=space_freed,
                    duration=0.3,
                    errors=[],
                )
            )
            total_items += items_deleted
            total_space += space_freed

        # Clean failed tasks
        if options.failed_tasks:
            items_deleted = 0

            if not options.dry_run:
                # Delete failed tasks older than cutoff
                deleted = (
                    self.db.query(PipelineTask)
                    .filter(and_(PipelineTask.status == "failed", PipelineTask.created_at < cutoff_date))
                    .delete()
                )
                self.db.commit()
                items_deleted = deleted

            results.append(
                CleanupResult(
                    operation="failed_tasks", items_deleted=items_deleted, space_freed=0, duration=0.8, errors=[]
                )
            )
            total_items += items_deleted

        # Clean old notifications
        if options.old_notifications:
            items_deleted = 0

            if not options.dry_run:
                # Delete read notifications older than cutoff
                deleted = (
                    self.db.query(Notification)
                    .filter(and_(Notification.read == True, Notification.created_at < cutoff_date))
                    .delete()
                )
                self.db.commit()
                items_deleted = deleted

            results.append(
                CleanupResult(
                    operation="old_notifications", items_deleted=items_deleted, space_freed=0, duration=0.4, errors=[]
                )
            )
            total_items += items_deleted

        duration = (datetime.utcnow() - start_time).total_seconds()

        return CleanupResponse(
            total_items_deleted=total_items,
            total_space_freed=total_space,
            total_duration=duration,
            results=results,
            timestamp=datetime.utcnow(),
        )

    # ========================================================================
    # System Statistics
    # ========================================================================

    def get_admin_system_stats(self) -> AdminSystemStats:
        """Get comprehensive system statistics"""

        # User statistics
        total_users = self.db.query(func.count(User.id)).scalar() or 0
        active_users = self.db.query(func.count(User.id)).filter(User.is_active == True).scalar() or 0
        admin_users = self.db.query(func.count(User.id)).filter(User.role == "admin").scalar() or 0

        # Resource statistics
        total_projects = self.db.query(func.count(Project.id)).scalar() or 0
        total_samples = self.db.query(func.count(Sample.id)).scalar() or 0
        total_tasks = self.db.query(func.count(PipelineTask.id)).scalar() or 0
        total_files = self.db.query(func.count(File.id)).scalar() or 0

        # Storage statistics
        total_storage_used = self.db.query(func.sum(User.storage_used)).scalar() or 0
        total_storage_allocated = self.db.query(func.sum(User.storage_quota)).scalar() or 0
        avg_storage_per_user = total_storage_used / total_users if total_users > 0 else 0

        # Task statistics
        today = datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0)
        week_ago = today - timedelta(days=7)
        month_ago = today - timedelta(days=30)

        tasks_today = self.db.query(func.count(PipelineTask.id)).filter(PipelineTask.created_at >= today).scalar() or 0

        tasks_this_week = (
            self.db.query(func.count(PipelineTask.id)).filter(PipelineTask.created_at >= week_ago).scalar() or 0
        )

        tasks_this_month = (
            self.db.query(func.count(PipelineTask.id)).filter(PipelineTask.created_at >= month_ago).scalar() or 0
        )

        # Success rates
        def calculate_success_rate(start_date):
            result = (
                self.db.query(
                    func.count(PipelineTask.id).label("total"),
                    func.sum(func.case((PipelineTask.status == "completed", 1), else_=0)).label("completed"),
                )
                .filter(PipelineTask.created_at >= start_date)
                .first()
            )

            if result and result.total > 0:
                return (result.completed / result.total) * 100
            return 0.0

        success_rate_today = calculate_success_rate(today)
        success_rate_this_week = calculate_success_rate(week_ago)
        success_rate_this_month = calculate_success_rate(month_ago)

        # System metrics (mock - in production would come from monitoring system)
        import psutil

        cpu_usage = psutil.cpu_percent(interval=1)
        memory = psutil.virtual_memory()
        disk = psutil.disk_usage("/")

        return AdminSystemStats(
            total_users=total_users,
            active_users=active_users,
            admin_users=admin_users,
            total_projects=total_projects,
            total_samples=total_samples,
            total_tasks=total_tasks,
            total_files=total_files,
            total_storage_used=total_storage_used,
            total_storage_allocated=total_storage_allocated,
            avg_storage_per_user=avg_storage_per_user,
            tasks_today=tasks_today,
            tasks_this_week=tasks_this_week,
            tasks_this_month=tasks_this_month,
            success_rate_today=success_rate_today,
            success_rate_this_week=success_rate_this_week,
            success_rate_this_month=success_rate_this_month,
            active_sessions=0,  # Would need session tracking
            cpu_usage=cpu_usage,
            memory_usage=memory.percent,
            disk_usage=disk.percent,
            generated_at=datetime.utcnow(),
        )

    # ========================================================================
    # Enhanced Features
    # ========================================================================

    def get_system_metrics(self) -> SystemMetrics:
        """Get detailed system metrics"""
        import psutil

        # CPU metrics
        cpu_percent = psutil.cpu_percent(interval=1)
        load_avg = psutil.getloadavg() if hasattr(psutil, "getloadavg") else [0, 0, 0]

        # Memory metrics
        memory = psutil.virtual_memory()

        # Disk metrics
        disk = psutil.disk_usage("/")

        # Network metrics (optional)
        try:
            network = psutil.net_io_counters()
            network_info = {
                "bytes_sent": network.bytes_sent,
                "bytes_recv": network.bytes_recv,
                "packets_sent": network.packets_sent,
                "packets_recv": network.packets_recv,
            }
        except Exception:
            network_info = None

        return SystemMetrics(
            cpu={"usage": cpu_percent, "load": list(load_avg)},
            memory={"used": memory.used, "total": memory.total, "usagePercent": memory.percent},
            disk={"used": disk.used, "total": disk.total, "usagePercent": disk.percent},
            network=network_info,
        )

    def get_alerts(self, resolved: bool = False) -> AlertListResponse:
        """Get system alerts from the database, generating new ones as needed."""
        from app.services.alert_service import AlertService

        alert_service = AlertService(self.db)

        # Run health check to generate any new alerts
        try:
            alert_service.check_system_health()
        except Exception as e:
            logger.error(f"Health check during alert query failed: {e}")

        # Query persisted alerts
        db_alerts = alert_service.get_alerts(resolved=resolved, limit=200)
        unresolved_count = alert_service.count_alerts(resolved=False)

        alerts = [
            Alert(
                id=str(a.id),
                type=AlertType(a.type),
                severity=AlertSeverity(a.severity),
                title=a.title,
                message=a.message,
                timestamp=a.timestamp,
                resolved=a.resolved,
                resolved_at=a.resolved_at,
                resolved_by=str(a.resolved_by) if a.resolved_by else None,
                source=a.source,
                metadata=a.alert_metadata,
            )
            for a in db_alerts
        ]

        return AlertListResponse(
            alerts=alerts,
            total=len(alerts),
            unresolved_count=unresolved_count,
        )

    def resolve_alert(self, alert_id: str, admin_user_id: str) -> Alert:
        """Resolve an alert in the database."""
        from app.services.alert_service import AlertService

        alert_service = AlertService(self.db)
        db_alert = alert_service.resolve_alert(alert_id, admin_user_id)

        if not db_alert:
            raise ValueError(f"Alert {alert_id} not found")

        return Alert(
            id=str(db_alert.id),
            type=AlertType(db_alert.type),
            severity=AlertSeverity(db_alert.severity),
            title=db_alert.title,
            message=db_alert.message,
            timestamp=db_alert.timestamp,
            resolved=db_alert.resolved,
            resolved_at=db_alert.resolved_at,
            resolved_by=str(db_alert.resolved_by) if db_alert.resolved_by else None,
            source=db_alert.source,
            metadata=db_alert.alert_metadata,
        )

    def get_audit_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        user_id: Optional[str] = None,
        action: Optional[str] = None,
        skip: int = 0,
        limit: int = 50,
    ) -> List[AuditLogEntry]:
        """Get audit logs from the database."""
        from app.services.audit_service import AuditService

        audit_service = AuditService(self.db)
        db_logs = audit_service.get_logs(
            start_date=start_date,
            end_date=end_date,
            admin_user_id=user_id,
            action=action,
            skip=skip,
            limit=limit,
        )

        return [
            AuditLogEntry(
                id=str(log.id),
                timestamp=log.timestamp,
                action=(
                    AuditAction(log.action)
                    if log.action in [a.value for a in AuditAction]
                    else AuditAction.CONFIG_UPDATED
                ),
                admin_user_id=str(log.admin_user_id),
                admin_username=log.admin_username,
                target_user_id=str(log.target_user_id) if log.target_user_id else None,
                target_username=log.target_username,
                details=log.details or {},
                ip_address=log.ip_address,
                user_agent=log.user_agent,
            )
            for log in db_logs
        ]

    def export_audit_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        user_id: Optional[str] = None,
        action: Optional[AuditAction] = None,
        format: str = "json",
    ) -> ExportResult:
        """Export audit logs to a downloadable file."""
        import json
        from pathlib import Path

        from app.services.audit_service import AuditService

        audit_service = AuditService(self.db)

        action_str = action.value if action else None
        logs = audit_service.get_logs(
            start_date=start_date,
            end_date=end_date,
            admin_user_id=user_id,
            action=action_str,
            skip=0,
            limit=100000,
        )

        # Determine export directory (use BACKUP_DIR/exports as a safe location)
        from app.core.config import settings

        export_dir = Path(settings.BACKUP_DIR) / "exports"
        export_dir.mkdir(parents=True, exist_ok=True)

        filename = f"audit_logs_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.{format}"
        file_path = export_dir / filename

        try:
            if format == "csv":
                csv_content = audit_service.export_logs_to_csv(
                    start_date=start_date,
                    end_date=end_date,
                    admin_user_id=user_id,
                    action=action_str,
                )
                file_path.write_text(csv_content)
            else:  # json
                data = [log.to_dict() for log in logs]
                file_path.write_text(json.dumps(data, indent=2, default=str))

            file_size = file_path.stat().st_size
        except Exception as e:
            logger.error(f"Failed to export audit logs: {e}")
            file_size = 0

        return ExportResult(
            download_url=f"/api/v1/admin/downloads/{filename}",
            file_name=filename,
            file_size=file_size,
            expires_at=datetime.utcnow() + timedelta(hours=24),
        )

    def get_resource_usage(self) -> ResourceUsage:
        """Get resource usage summary using real metrics."""
        import psutil

        from app.services.job_service import JobService

        # Storage usage
        disk = psutil.disk_usage("/")
        storage_used = self.db.query(func.sum(User.storage_used)).scalar() or 0

        # Compute usage from actual job queue
        job_service = JobService(self.db)
        active_jobs = job_service.get_active_jobs_count()

        # Active pipeline tasks
        active_tasks = (
            self.db.query(func.count(PipelineTask.id)).filter(PipelineTask.status.in_(["pending", "running"])).scalar()
            or 0
        )

        # Compute limit from CPU count (rough estimate)
        cpu_count = psutil.cpu_count() or 4
        compute_limit = cpu_count * 4  # Allow 4 jobs per CPU

        # Memory usage
        memory = psutil.virtual_memory()

        return ResourceUsage(
            storage={"total": int(disk.total), "used": int(storage_used)},
            compute={"active": int(active_jobs + active_tasks), "limit": compute_limit},
            memory={"total": int(memory.total), "used": int(memory.used)},
        )

    def create_backup(
        self, backup_type: BackupType, description: Optional[str], admin_user_id: str, compress: bool = True
    ) -> BackupInfo:
        """
        Create system backup asynchronously via Celery.

        Returns a SystemBackup record with status='pending' that will be
        updated by the Celery task as the backup progresses.
        """
        from app.core.config import settings as app_settings
        from app.models.backup import SystemBackup
        from app.services.job_service import JobService

        # Create initial backup record
        backup = SystemBackup(
            backup_type=backup_type.value,
            status="pending",
            compressed=compress,
            description=description,
            created_by=admin_user_id,
            expires_at=datetime.utcnow() + timedelta(days=app_settings.BACKUP_RETENTION_DAYS),
        )
        self.db.add(backup)
        self.db.commit()
        self.db.refresh(backup)

        # Dispatch Celery task and create corresponding job record
        try:
            from app.workers.admin_tasks import create_backup_task

            celery_result = create_backup_task.delay(
                backup_id=str(backup.id),
                backup_type=backup_type.value,
                admin_user_id=admin_user_id,
                description=description,
                compress=compress,
            )

            # Track job
            job_service = JobService(self.db)
            job_service.create_job(
                type="backup",
                user_id=admin_user_id,
                celery_task_id=celery_result.id,
                parameters={
                    "backup_id": str(backup.id),
                    "backup_type": backup_type.value,
                    "compress": compress,
                },
                message=f"Creating {backup_type.value} backup",
            )
        except Exception as e:
            logger.error(f"Failed to dispatch backup task: {e}")
            # Fallback: run synchronously
            from app.services.backup_service import BackupService

            backup_service = BackupService(self.db)
            backup = backup_service.create_backup(
                backup_type=backup_type.value,
                admin_user_id=admin_user_id,
                description=description,
                compress=compress,
            )

        return BackupInfo(
            id=str(backup.id),
            backup_type=BackupType(backup.backup_type),
            file_path=backup.file_path or "",
            size=backup.size or 0,
            compressed=backup.compressed,
            description=backup.description,
            created_at=backup.created_at,
            created_by=str(backup.created_by),
            status=backup.status,
        )

    def list_backups(self) -> BackupListResponse:
        """List all backups from the database."""
        from app.services.backup_service import BackupService

        backup_service = BackupService(self.db)
        db_backups = backup_service.list_backups(limit=100)

        backups = [
            BackupInfo(
                id=str(b.id),
                backup_type=BackupType(b.backup_type),
                file_path=b.file_path or "",
                size=b.size or 0,
                compressed=b.compressed,
                description=b.description,
                created_at=b.created_at,
                created_by=str(b.created_by),
                status=b.status,
            )
            for b in db_backups
        ]

        return BackupListResponse(
            backups=backups,
            total=len(backups),
        )

    def list_jobs(self, status: Optional[JobStatus] = None, skip: int = 0, limit: int = 50) -> JobListResponse:
        """List system jobs from the database, syncing with Celery state."""
        from app.services.job_service import JobService

        job_service = JobService(self.db)
        status_str = status.value if status else None
        db_jobs = job_service.list_jobs(status=status_str, skip=skip, limit=limit)

        # Sync running jobs with Celery to ensure up-to-date status
        for job in db_jobs:
            if job.status in ("pending", "running") and job.celery_task_id:
                try:
                    job_service.sync_with_celery(job)
                except Exception as e:
                    logger.error(f"Failed to sync job {job.id} with Celery: {e}")

        total = job_service.count_jobs(status=status_str)

        jobs = [
            JobInfo(
                id=str(j.id),
                type=JobType(j.type) if j.type in [t.value for t in JobType] else JobType.SYSTEM_TASK,
                status=JobStatus(j.status),
                user_id=str(j.user_id) if j.user_id else "",
                username=j.username,
                created_at=j.created_at,
                started_at=j.started_at,
                completed_at=j.completed_at,
                progress=j.progress,
                message=j.message,
                result=j.result,
                error=j.error,
            )
            for j in db_jobs
        ]

        return JobListResponse(
            jobs=jobs,
            total=total,
            page=skip // limit + 1 if limit > 0 else 1,
            page_size=limit,
        )

    def cancel_job(self, job_id: str) -> JobInfo:
        """Cancel a running job via Celery and update database state."""
        from app.services.job_service import JobService

        job_service = JobService(self.db)
        job = job_service.cancel_job(job_id)

        if not job:
            raise ValueError(f"Job {job_id} not found")

        return JobInfo(
            id=str(job.id),
            type=JobType(job.type) if job.type in [t.value for t in JobType] else JobType.SYSTEM_TASK,
            status=JobStatus(job.status),
            user_id=str(job.user_id) if job.user_id else "",
            username=job.username,
            created_at=job.created_at,
            started_at=job.started_at,
            completed_at=job.completed_at,
            progress=job.progress,
            message=job.message,
            result=job.result,
            error=job.error,
        )

    def retry_job(self, job_id: str) -> JobInfo:
        """Retry a failed job by creating a new job record and dispatching it."""
        from app.services.job_service import JobService

        job_service = JobService(self.db)
        new_job = job_service.retry_job(job_id)

        if not new_job:
            raise ValueError(f"Job {job_id} not found or not in failed state")

        # Re-dispatch the job to Celery based on its type
        try:
            if new_job.type == "backup":
                from app.workers.admin_tasks import create_backup_task

                params = new_job.parameters or {}
                celery_result = create_backup_task.delay(
                    backup_id=params.get("backup_id"),
                    backup_type=params.get("backup_type", "full"),
                    admin_user_id=str(new_job.user_id) if new_job.user_id else None,
                    description=params.get("description"),
                    compress=params.get("compress", True),
                )
                new_job.celery_task_id = celery_result.id
                self.db.commit()
                self.db.refresh(new_job)
        except Exception as e:
            logger.error(f"Failed to dispatch retry for job {new_job.id}: {e}")

        return JobInfo(
            id=str(new_job.id),
            type=JobType(new_job.type) if new_job.type in [t.value for t in JobType] else JobType.SYSTEM_TASK,
            status=JobStatus(new_job.status),
            user_id=str(new_job.user_id) if new_job.user_id else "",
            username=new_job.username,
            created_at=new_job.created_at,
            started_at=new_job.started_at,
            completed_at=new_job.completed_at,
            progress=new_job.progress,
            message=new_job.message,
        )
