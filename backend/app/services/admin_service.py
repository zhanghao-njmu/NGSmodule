"""
Admin service for user management, system configuration, and logs
"""
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_, desc
from datetime import datetime, timedelta
from typing import List, Dict, Optional, Any, Tuple
import uuid
import os
import json
import logging
from pathlib import Path
import math

from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.models.file import File
from app.models.notification import Notification
from app.schemas.admin import *
from app.core.security import get_password_hash


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
        sort_order: str = "desc"
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
                    User.organization.ilike(search_term)
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
                    updated_at=user.updated_at
                )
            )

        total_pages = math.ceil(total / limit) if limit > 0 else 0
        page = (skip // limit) + 1 if limit > 0 else 1

        return UserListResponse(
            users=user_list,
            total=total,
            page=page,
            page_size=limit,
            total_pages=total_pages
        )

    def get_user_detail(self, user_id: str) -> Optional[AdminUserDetail]:
        """Get detailed information about a user"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return None

        # Count user's resources
        total_projects = self.db.query(func.count(Project.id)).filter(
            Project.owner_id == user_id
        ).scalar() or 0

        total_samples = self.db.query(func.count(Sample.id)).join(Project).filter(
            Project.owner_id == user_id
        ).scalar() or 0

        total_tasks = self.db.query(func.count(PipelineTask.id)).join(Sample).join(Project).filter(
            Project.owner_id == user_id
        ).scalar() or 0

        total_files = self.db.query(func.count(File.id)).join(Sample).join(Project).filter(
            Project.owner_id == user_id
        ).scalar() or 0

        # Get last activity
        last_activity = self.db.query(func.max(PipelineTask.updated_at)).join(Sample).join(Project).filter(
            Project.owner_id == user_id
        ).scalar()

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
            login_count=0  # Would need login tracking table
        )

    def update_user(
        self,
        user_id: str,
        update_data: UserUpdateRequest
    ) -> Optional[AdminUserDetail]:
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
        self,
        user_id: str,
        is_active: bool,
        reason: Optional[str] = None
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

    def reset_user_password(
        self,
        user_id: str,
        new_password: str,
        notify_user: bool = True
    ) -> bool:
        """Reset user password"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False

        user.password_hash = get_password_hash(new_password)
        user.updated_at = datetime.utcnow()

        self.db.commit()

        logger.info(f"Password reset for user {user.username}")

        # TODO: Send email notification if notify_user is True

        return True

    def delete_user(
        self,
        user_id: str,
        transfer_data_to: Optional[str] = None
    ) -> bool:
        """Delete a user (with optional data transfer)"""

        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            return False

        # Transfer data if requested
        if transfer_data_to:
            target_user = self.db.query(User).filter(User.id == transfer_data_to).first()
            if target_user:
                # Transfer projects
                self.db.query(Project).filter(Project.owner_id == user_id).update(
                    {"owner_id": transfer_data_to}
                )
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
            last_updated=datetime.utcnow()
        )

    def update_system_config(
        self,
        category: SystemConfigCategory,
        updates: Dict[str, Any],
        admin_id: str
    ) -> SystemConfig:
        """Update system configuration"""

        # In production, this would update a config table or file
        logger.info(f"System config {category.value} updated by admin {admin_id}: {updates}")

        # TODO: Implement actual config persistence

        return self.get_system_config()

    def reset_system_config(
        self,
        categories: Optional[List[SystemConfigCategory]] = None
    ) -> SystemConfig:
        """Reset system configuration to defaults"""

        if categories:
            logger.info(f"Resetting config categories: {[c.value for c in categories]}")
        else:
            logger.info("Resetting all config to defaults")

        # TODO: Implement actual config reset

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
        offset: int = 0
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
                    ip_address="192.168.1.1"
                )
            )

        # Apply search filter
        if search:
            logs = [log for log in logs if search.lower() in log.message.lower()]

        total = len(logs)
        has_more = total > (offset + limit)

        # Apply pagination
        logs = logs[offset:offset + limit]

        return LogResponse(
            logs=logs,
            total=total,
            has_more=has_more
        )

    def download_logs(
        self,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        levels: Optional[List[LogLevel]] = None,
        sources: Optional[List[LogSource]] = None,
        format: str = "json"
    ) -> str:
        """Generate log file for download"""

        logs_response = self.get_logs(
            start_date=start_date,
            end_date=end_date,
            levels=levels,
            sources=sources,
            limit=10000
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
                [f"[{log.timestamp}] {log.level.value.upper()} [{log.source.value}] {log.message}"
                 for log in logs_response.logs]
            )

        # Save to temp file and return path
        temp_file = f"/tmp/logs_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.{format}"
        with open(temp_file, 'w') as f:
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
                    message="Database is responsive"
                )
            )
        except Exception as e:
            services.append(
                ServiceHealth(
                    name="PostgreSQL",
                    status=ServiceStatus.DOWN,
                    last_check=datetime.utcnow(),
                    message=f"Database error: {str(e)}"
                )
            )

        # Check Redis (mock)
        services.append(
            ServiceHealth(
                name="Redis",
                status=ServiceStatus.HEALTHY,
                response_time=1.5,
                last_check=datetime.utcnow(),
                message="Cache is operational"
            )
        )

        # Check MinIO (mock)
        services.append(
            ServiceHealth(
                name="MinIO",
                status=ServiceStatus.HEALTHY,
                response_time=8.3,
                last_check=datetime.utcnow(),
                message="Object storage is available"
            )
        )

        # Check Celery (mock)
        services.append(
            ServiceHealth(
                name="Celery",
                status=ServiceStatus.HEALTHY,
                response_time=3.1,
                last_check=datetime.utcnow(),
                message="Task queue is processing"
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
            version="1.0.0"
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
                    operation="old_logs",
                    items_deleted=items_deleted,
                    space_freed=space_freed,
                    duration=0.5,
                    errors=[]
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
                    errors=[]
                )
            )
            total_items += items_deleted
            total_space += space_freed

        # Clean failed tasks
        if options.failed_tasks:
            items_deleted = 0

            if not options.dry_run:
                # Delete failed tasks older than cutoff
                deleted = self.db.query(PipelineTask).filter(
                    and_(
                        PipelineTask.status == "failed",
                        PipelineTask.created_at < cutoff_date
                    )
                ).delete()
                self.db.commit()
                items_deleted = deleted

            results.append(
                CleanupResult(
                    operation="failed_tasks",
                    items_deleted=items_deleted,
                    space_freed=0,
                    duration=0.8,
                    errors=[]
                )
            )
            total_items += items_deleted

        # Clean old notifications
        if options.old_notifications:
            items_deleted = 0

            if not options.dry_run:
                # Delete read notifications older than cutoff
                deleted = self.db.query(Notification).filter(
                    and_(
                        Notification.read == True,
                        Notification.created_at < cutoff_date
                    )
                ).delete()
                self.db.commit()
                items_deleted = deleted

            results.append(
                CleanupResult(
                    operation="old_notifications",
                    items_deleted=items_deleted,
                    space_freed=0,
                    duration=0.4,
                    errors=[]
                )
            )
            total_items += items_deleted

        duration = (datetime.utcnow() - start_time).total_seconds()

        return CleanupResponse(
            total_items_deleted=total_items,
            total_space_freed=total_space,
            total_duration=duration,
            results=results,
            timestamp=datetime.utcnow()
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

        tasks_today = self.db.query(func.count(PipelineTask.id)).filter(
            PipelineTask.created_at >= today
        ).scalar() or 0

        tasks_this_week = self.db.query(func.count(PipelineTask.id)).filter(
            PipelineTask.created_at >= week_ago
        ).scalar() or 0

        tasks_this_month = self.db.query(func.count(PipelineTask.id)).filter(
            PipelineTask.created_at >= month_ago
        ).scalar() or 0

        # Success rates
        def calculate_success_rate(start_date):
            result = self.db.query(
                func.count(PipelineTask.id).label('total'),
                func.sum(func.case((PipelineTask.status == 'completed', 1), else_=0)).label('completed')
            ).filter(PipelineTask.created_at >= start_date).first()

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
        disk = psutil.disk_usage('/')

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
            generated_at=datetime.utcnow()
        )
