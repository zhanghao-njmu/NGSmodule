"""
Statistics service for aggregating and calculating system statistics
"""

from datetime import datetime, timedelta
from typing import Optional

import psutil
from sqlalchemy import func
from sqlalchemy.orm import Session

from app.models.file import File
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask as Task
from app.models.user import User
from app.schemas.stats import (
    FileStats,
    PipelineStats,
    ProjectStats,
    QuickStats,
    SampleStats,
    StatsSummary,
    StorageStats,
    SystemStats,
    TaskStats,
    TrendData,
    UserActivityStats,
)


class StatsService:
    """Service for calculating and retrieving statistics"""

    def __init__(self, db: Session):
        self.db = db

    def get_summary(self, include_system: bool = False) -> StatsSummary:
        """
        Get comprehensive statistics summary

        Args:
            include_system: Whether to include system-level stats (admin only)

        Returns:
            StatsSummary with all statistics
        """
        return StatsSummary(
            projects=self.get_project_stats(),
            samples=self.get_sample_stats(),
            tasks=self.get_task_stats(),
            files=self.get_file_stats(),
            storage=self.get_storage_stats(),
            users=self.get_user_activity_stats() if include_system else None,
            pipelines=self.get_pipeline_stats() if include_system else None,
            system=self.get_system_stats() if include_system else None,
        )

    def get_project_stats(self) -> ProjectStats:
        """Get project statistics"""
        total = self.db.query(func.count(Project.id)).scalar() or 0

        active = self.db.query(func.count(Project.id)).filter(Project.status.in_(["active", "running"])).scalar() or 0

        completed = self.db.query(func.count(Project.id)).filter(Project.status == "completed").scalar() or 0

        failed = self.db.query(func.count(Project.id)).filter(Project.status == "failed").scalar() or 0

        return ProjectStats(total=total, active=active, completed=completed, failed=failed)

    def get_sample_stats(self) -> SampleStats:
        """Get sample statistics"""
        total = self.db.query(func.count(Sample.id)).scalar() or 0

        processing = (
            self.db.query(func.count(Sample.id)).filter(Sample.status.in_(["processing", "running"])).scalar() or 0
        )

        completed = self.db.query(func.count(Sample.id)).filter(Sample.status == "completed").scalar() or 0

        failed = self.db.query(func.count(Sample.id)).filter(Sample.status == "failed").scalar() or 0

        return SampleStats(total=total, processing=processing, completed=completed, failed=failed)

    def get_task_stats(self) -> TaskStats:
        """Get task statistics"""
        total = self.db.query(func.count(Task.id)).scalar() or 0

        pending = self.db.query(func.count(Task.id)).filter(Task.status == "pending").scalar() or 0

        running = self.db.query(func.count(Task.id)).filter(Task.status == "running").scalar() or 0

        completed = self.db.query(func.count(Task.id)).filter(Task.status == "completed").scalar() or 0

        failed = self.db.query(func.count(Task.id)).filter(Task.status == "failed").scalar() or 0

        cancelled = self.db.query(func.count(Task.id)).filter(Task.status == "cancelled").scalar() or 0

        return TaskStats(
            total=total, pending=pending, running=running, completed=completed, failed=failed, cancelled=cancelled
        )

    def get_file_stats(self) -> FileStats:
        """Get file statistics"""
        total = self.db.query(func.count(File.id)).scalar() or 0
        total_size = self.db.query(func.sum(File.size)).scalar() or 0

        # Group files by type
        file_types = self.db.query(File.file_type, func.count(File.id)).group_by(File.file_type).all()

        by_type = {file_type: count for file_type, count in file_types}

        return FileStats(total=total, total_size=int(total_size), by_type=by_type)

    def get_storage_stats(self) -> StorageStats:
        """Get storage usage statistics"""
        # Get total file size from database
        total_files_size = self.db.query(func.sum(File.size)).scalar() or 0
        files_count = self.db.query(func.count(File.id)).scalar() or 0

        # Try to get actual disk usage
        try:
            disk = psutil.disk_usage("/")
            total_space = disk.total
            used_space = disk.used
            available_space = disk.free
            usage_percentage = (used_space / total_space) * 100
        except Exception:
            # Fallback if psutil not available
            total_space = 1000000000000  # 1TB default
            used_space = int(total_files_size)
            available_space = total_space - used_space
            usage_percentage = (used_space / total_space) * 100

        # Calculate space usage by category
        projects_space = (
            self.db.query(func.sum(File.size))
            .join(Sample, File.sample_id == Sample.id)
            .join(Project, Sample.project_id == Project.id)
            .scalar()
            or 0
        )

        return StorageStats(
            total_space=total_space,
            used_space=used_space,
            available_space=available_space,
            usage_percentage=round(usage_percentage, 2),
            files_count=files_count,
            projects_space=int(projects_space),
            samples_space=int(total_files_size),
            results_space=0,  # Can be calculated if needed
        )

    def get_user_activity_stats(self) -> UserActivityStats:
        """Get user activity statistics"""
        now = datetime.utcnow()
        today = now.replace(hour=0, minute=0, second=0, microsecond=0)
        week_ago = now - timedelta(days=7)
        month_ago = now - timedelta(days=30)

        total_users = self.db.query(func.count(User.id)).scalar() or 0

        # Active users (users with recent login)
        active_today = self.db.query(func.count(User.id)).filter(User.last_login >= today).scalar() or 0

        active_week = self.db.query(func.count(User.id)).filter(User.last_login >= week_ago).scalar() or 0

        active_month = self.db.query(func.count(User.id)).filter(User.last_login >= month_ago).scalar() or 0

        # New users
        new_week = self.db.query(func.count(User.id)).filter(User.created_at >= week_ago).scalar() or 0

        new_month = self.db.query(func.count(User.id)).filter(User.created_at >= month_ago).scalar() or 0

        return UserActivityStats(
            total_users=total_users,
            active_users_today=active_today,
            active_users_week=active_week,
            active_users_month=active_month,
            new_users_week=new_week,
            new_users_month=new_month,
        )

    def get_pipeline_stats(self) -> PipelineStats:
        """Get pipeline execution statistics"""
        total_runs = self.db.query(func.count(Task.id)).scalar() or 0

        completed_runs = self.db.query(func.count(Task.id)).filter(Task.status == "completed").scalar() or 0

        success_rate = (completed_runs / total_runs * 100) if total_runs > 0 else 0

        # Calculate average duration for completed tasks
        avg_duration = (
            self.db.query(func.avg(func.extract("epoch", Task.completed_at - Task.started_at)))
            .filter(Task.status == "completed", Task.started_at.isnot(None), Task.completed_at.isnot(None))
            .scalar()
            or 0
        )

        # Most used pipelines
        most_used = (
            self.db.query(Task.pipeline_id, func.count(Task.id).label("count"))
            .filter(Task.pipeline_id.isnot(None))
            .group_by(Task.pipeline_id)
            .order_by(func.count(Task.id).desc())
            .limit(5)
            .all()
        )

        most_used_list = [{"pipeline_id": str(pipeline_id), "count": count} for pipeline_id, count in most_used]

        return PipelineStats(
            total_runs=total_runs,
            success_rate=round(success_rate, 2),
            average_duration=round(float(avg_duration), 2),
            most_used=most_used_list,
        )

    def get_system_stats(self) -> SystemStats:
        """Get system-level statistics"""
        # These would typically come from monitoring system
        # For now, return placeholder values

        try:
            # Get system uptime
            boot_time = psutil.boot_time()
            uptime = datetime.now().timestamp() - boot_time
        except Exception:
            uptime = 0

        return SystemStats(
            uptime=uptime,
            api_calls_today=0,  # Requires logging/monitoring system
            api_calls_hour=0,  # Requires logging/monitoring system
            active_connections=0,  # Requires WebSocket tracking
            queue_size=0,  # Requires Celery queue inspection
        )

    def get_quick_stats(self, user_id: Optional[str] = None) -> QuickStats:
        """
        Get quick statistics for header/widget display

        Args:
            user_id: Optional user ID to filter user-specific stats

        Returns:
            QuickStats with essential metrics
        """
        query_filter = []
        if user_id:
            query_filter.append(Project.user_id == user_id)

        total_projects = self.db.query(func.count(Project.id)).filter(*query_filter).scalar() or 0

        total_samples = (
            self.db.query(func.count(Sample.id))
            .join(Project, Sample.project_id == Project.id)
            .filter(*query_filter)
            .scalar()
            or 0
        )

        running_tasks = self.db.query(func.count(Task.id)).filter(Task.status == "running").scalar() or 0

        # Storage stats
        storage_used = self.db.query(func.sum(File.size)).scalar() or 0
        storage_quota = 1000000000000  # 1TB default, should come from user settings

        return QuickStats(
            total_projects=total_projects,
            total_samples=total_samples,
            running_tasks=running_tasks,
            storage_used=int(storage_used),
            storage_quota=storage_quota,
        )

    def get_trend_data(self, metric: str, period: str = "week", days: int = 7) -> TrendData:
        """
        Get trend data for time-series analysis

        Args:
            metric: Metric name (projects, samples, tasks, etc.)
            period: Time period (day, week, month)
            days: Number of days to include

        Returns:
            TrendData with time-series data points
        """
        now = datetime.utcnow()
        start_date = now - timedelta(days=days)

        data_points = []

        for i in range(days):
            day = start_date + timedelta(days=i)
            next_day = day + timedelta(days=1)

            if metric == "projects":
                count = (
                    self.db.query(func.count(Project.id))
                    .filter(Project.created_at >= day, Project.created_at < next_day)
                    .scalar()
                    or 0
                )
            elif metric == "samples":
                count = (
                    self.db.query(func.count(Sample.id))
                    .filter(Sample.created_at >= day, Sample.created_at < next_day)
                    .scalar()
                    or 0
                )
            elif metric == "tasks":
                count = (
                    self.db.query(func.count(Task.id))
                    .filter(Task.created_at >= day, Task.created_at < next_day)
                    .scalar()
                    or 0
                )
            else:
                count = 0

            data_points.append({"date": day.isoformat(), "value": count})

        return TrendData(metric=metric, period=period, data_points=data_points)
