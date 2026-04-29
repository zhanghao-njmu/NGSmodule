"""
Analytics service for data analysis and visualization
"""

from datetime import datetime, timedelta
from statistics import mean, stdev
from typing import List, Optional

from sqlalchemy import and_, case, desc, func
from sqlalchemy.orm import Session

from app.models.file import File
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.models.user import User
from app.schemas.analytics import *


class AnalyticsService:
    """Service for analytics operations"""

    def __init__(self, db: Session):
        self.db = db

    # ========================================================================
    # Time Series Analytics
    # ========================================================================

    def get_time_series_data(
        self,
        metric: str,
        time_range: TimeRange,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        user_id: Optional[str] = None,
    ) -> TimeSeriesData:
        """Get time series data for a metric"""

        # Default time range
        if not end_date:
            end_date = datetime.utcnow()
        if not start_date:
            if time_range == TimeRange.DAY:
                start_date = end_date - timedelta(days=1)
            elif time_range == TimeRange.WEEK:
                start_date = end_date - timedelta(weeks=1)
            elif time_range == TimeRange.MONTH:
                start_date = end_date - timedelta(days=30)
            elif time_range == TimeRange.YEAR:
                start_date = end_date - timedelta(days=365)
            else:
                start_date = end_date - timedelta(days=7)

        # Get data based on metric
        if metric == "tasks":
            data_points = self._get_task_time_series(start_date, end_date, user_id)
        elif metric == "samples":
            data_points = self._get_sample_time_series(start_date, end_date, user_id)
        elif metric == "storage":
            data_points = self._get_storage_time_series(start_date, end_date, user_id)
        else:
            data_points = []

        # Calculate statistics
        values = [dp.value for dp in data_points]
        statistics = {}
        if values:
            statistics = {
                "min": min(values),
                "max": max(values),
                "avg": mean(values),
                "std": stdev(values) if len(values) > 1 else 0.0,
            }

        return TimeSeriesData(
            metric=metric,
            time_range=time_range,
            data_points=data_points,
            total_points=len(data_points),
            start_time=start_date,
            end_time=end_date,
            statistics=statistics,
        )

    def _get_task_time_series(
        self, start_date: datetime, end_date: datetime, user_id: Optional[str]
    ) -> List[TimeSeriesDataPoint]:
        """Get task count time series"""

        query = self.db.query(
            func.date_trunc("day", PipelineTask.created_at).label("date"), func.count(PipelineTask.id).label("count")
        ).filter(and_(PipelineTask.created_at >= start_date, PipelineTask.created_at <= end_date))

        if user_id:
            # Filter by user's projects
            query = query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        results = query.group_by("date").order_by("date").all()

        return [
            TimeSeriesDataPoint(timestamp=row.date, value=float(row.count), label=row.date.strftime("%Y-%m-%d"))
            for row in results
        ]

    def _get_sample_time_series(
        self, start_date: datetime, end_date: datetime, user_id: Optional[str]
    ) -> List[TimeSeriesDataPoint]:
        """Get sample count time series"""

        query = self.db.query(
            func.date_trunc("day", Sample.created_at).label("date"), func.count(Sample.id).label("count")
        ).filter(and_(Sample.created_at >= start_date, Sample.created_at <= end_date))

        if user_id:
            query = query.join(Project).filter(Project.owner_id == user_id)

        results = query.group_by("date").order_by("date").all()

        return [
            TimeSeriesDataPoint(timestamp=row.date, value=float(row.count), label=row.date.strftime("%Y-%m-%d"))
            for row in results
        ]

    def _get_storage_time_series(
        self, start_date: datetime, end_date: datetime, user_id: Optional[str]
    ) -> List[TimeSeriesDataPoint]:
        """Get storage usage time series"""

        query = self.db.query(
            func.date_trunc("day", File.created_at).label("date"), func.sum(File.file_size).label("total_size")
        ).filter(and_(File.created_at >= start_date, File.created_at <= end_date))

        if user_id:
            query = query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        results = query.group_by("date").order_by("date").all()

        # Calculate cumulative storage
        cumulative = 0
        data_points = []
        for row in results:
            cumulative += row.total_size or 0
            data_points.append(
                TimeSeriesDataPoint(timestamp=row.date, value=float(cumulative), label=row.date.strftime("%Y-%m-%d"))
            )

        return data_points

    # ========================================================================
    # Project Analytics
    # ========================================================================

    def get_project_performance(
        self, project_id: str, user_id: Optional[str] = None
    ) -> Optional[ProjectPerformanceMetrics]:
        """Get performance metrics for a project"""

        project = self.db.query(Project).filter(Project.id == project_id).first()
        if not project:
            return None

        # Check permission
        if user_id and project.owner_id != user_id:
            return None

        # Count samples
        total_samples = self.db.query(func.count(Sample.id)).filter(Sample.project_id == project_id).scalar()

        # Count tasks
        task_stats = (
            self.db.query(
                func.count(PipelineTask.id).label("total"),
                func.sum(case((PipelineTask.status == "completed", 1), else_=0)).label("completed"),
                func.sum(case((PipelineTask.status == "failed", 1), else_=0)).label("failed"),
                func.avg(func.extract("epoch", PipelineTask.end_time - PipelineTask.start_time)).label("avg_duration"),
            )
            .join(Sample)
            .filter(
                Sample.project_id == project_id, PipelineTask.start_time.isnot(None), PipelineTask.end_time.isnot(None)
            )
            .first()
        )

        total_tasks = task_stats.total or 0
        completed_tasks = task_stats.completed or 0
        failed_tasks = task_stats.failed or 0
        success_rate = (completed_tasks / total_tasks * 100) if total_tasks > 0 else 0.0

        # Calculate total processing time
        total_processing = (
            self.db.query(func.sum(func.extract("epoch", PipelineTask.end_time - PipelineTask.start_time)))
            .join(Sample)
            .filter(
                Sample.project_id == project_id, PipelineTask.start_time.isnot(None), PipelineTask.end_time.isnot(None)
            )
            .scalar()
        )

        # Calculate storage used
        storage_used = (
            self.db.query(func.sum(File.file_size)).join(Sample).filter(Sample.project_id == project_id).scalar() or 0
        )

        # Get last activity
        last_activity = (
            self.db.query(func.max(PipelineTask.updated_at))
            .join(Sample)
            .filter(Sample.project_id == project_id)
            .scalar()
        )

        return ProjectPerformanceMetrics(
            project_id=str(project.id),
            project_name=project.name,
            total_samples=total_samples,
            total_tasks=total_tasks,
            completed_tasks=completed_tasks,
            failed_tasks=failed_tasks,
            success_rate=success_rate,
            avg_task_duration=task_stats.avg_duration,
            total_processing_time=total_processing,
            storage_used=storage_used,
            created_at=project.created_at,
            last_activity=last_activity,
        )

    def compare_projects(self, project_ids: List[str], metric: str, user_id: Optional[str] = None) -> ProjectComparison:
        """Compare multiple projects"""

        projects = []
        for project_id in project_ids:
            metrics = self.get_project_performance(project_id, user_id)
            if metrics:
                projects.append(metrics)

        return ProjectComparison(projects=projects, comparison_metric=metric, generated_at=datetime.utcnow())

    # ========================================================================
    # Sample Quality Analytics
    # ========================================================================

    def get_sample_quality_distribution(
        self, metric_name: str, project_id: Optional[str] = None, user_id: Optional[str] = None
    ) -> SampleQualityDistribution:
        """Get distribution of sample quality metrics"""

        # This is a placeholder - in production, quality metrics would come from
        # pipeline results or QC reports

        query = self.db.query(Sample)

        if project_id:
            query = query.filter(Sample.project_id == project_id)

        if user_id:
            query = query.join(Project).filter(Project.owner_id == user_id)

        samples = query.all()
        total_samples = len(samples)

        # Mock distribution bins
        bins = [
            {"range": "0-20", "count": int(total_samples * 0.1)},
            {"range": "20-40", "count": int(total_samples * 0.15)},
            {"range": "40-60", "count": int(total_samples * 0.25)},
            {"range": "60-80", "count": int(total_samples * 0.35)},
            {"range": "80-100", "count": int(total_samples * 0.15)},
        ]

        statistics = {"min": 0.0, "max": 100.0, "avg": 65.5, "median": 70.0, "std": 15.2}

        return SampleQualityDistribution(
            metric_name=metric_name, bins=bins, total_samples=total_samples, statistics=statistics
        )

    # ========================================================================
    # Task Execution Analytics
    # ========================================================================

    def get_pipeline_performance(
        self, pipeline_name: Optional[str] = None, user_id: Optional[str] = None
    ) -> List[PipelinePerformanceMetrics]:
        """Get pipeline performance metrics"""

        query = self.db.query(
            PipelineTask.pipeline,
            func.count(PipelineTask.id).label("total"),
            func.sum(case((PipelineTask.status == "completed", 1), else_=0)).label("successful"),
            func.sum(case((PipelineTask.status == "failed", 1), else_=0)).label("failed"),
            func.avg(func.extract("epoch", PipelineTask.end_time - PipelineTask.start_time)).label("avg_duration"),
            func.min(func.extract("epoch", PipelineTask.end_time - PipelineTask.start_time)).label("min_duration"),
            func.max(func.extract("epoch", PipelineTask.end_time - PipelineTask.start_time)).label("max_duration"),
        ).filter(PipelineTask.start_time.isnot(None), PipelineTask.end_time.isnot(None))

        if pipeline_name:
            query = query.filter(PipelineTask.pipeline == pipeline_name)

        if user_id:
            query = query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        results = query.group_by(PipelineTask.pipeline).all()

        metrics = []
        for row in results:
            total = row.total or 0
            successful = row.successful or 0
            success_rate = (successful / total * 100) if total > 0 else 0.0

            metrics.append(
                PipelinePerformanceMetrics(
                    pipeline_name=row.pipeline,
                    total_runs=total,
                    successful_runs=successful,
                    failed_runs=row.failed or 0,
                    success_rate=success_rate,
                    avg_duration=row.avg_duration,
                    min_duration=row.min_duration,
                    max_duration=row.max_duration,
                    avg_cpu_usage=None,  # Would come from monitoring system
                    avg_memory_usage=None,
                )
            )

        return metrics

    def get_task_execution_trend(
        self,
        time_range: TimeRange,
        start_date: Optional[datetime] = None,
        end_date: Optional[datetime] = None,
        user_id: Optional[str] = None,
    ) -> TaskExecutionTrend:
        """Get task execution trend over time"""

        if not end_date:
            end_date = datetime.utcnow()
        if not start_date:
            if time_range == TimeRange.WEEK:
                start_date = end_date - timedelta(weeks=1)
            elif time_range == TimeRange.MONTH:
                start_date = end_date - timedelta(days=30)
            else:
                start_date = end_date - timedelta(days=7)

        query = self.db.query(
            func.date_trunc("day", PipelineTask.created_at).label("date"),
            func.count(PipelineTask.id).label("total"),
            func.sum(case((PipelineTask.status == "completed", 1), else_=0)).label("completed"),
            func.sum(case((PipelineTask.status == "failed", 1), else_=0)).label("failed"),
            func.sum(case((PipelineTask.status == "running", 1), else_=0)).label("running"),
        ).filter(and_(PipelineTask.created_at >= start_date, PipelineTask.created_at <= end_date))

        if user_id:
            query = query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        results = query.group_by("date").order_by("date").all()

        data = [
            {
                "date": row.date.strftime("%Y-%m-%d"),
                "total": row.total,
                "completed": row.completed or 0,
                "failed": row.failed or 0,
                "running": row.running or 0,
            }
            for row in results
        ]

        return TaskExecutionTrend(period=time_range, data=data, total_periods=len(data), generated_at=datetime.utcnow())

    # ========================================================================
    # Resource Analytics
    # ========================================================================

    def get_storage_analytics(self, user_id: Optional[str] = None) -> StorageAnalytics:
        """Get storage analytics"""

        if user_id:
            # Get user's storage quota and usage
            user = self.db.query(User).filter(User.id == user_id).first()
            if not user:
                return None

            total_storage = user.storage_quota
            used_storage = user.storage_used
        else:
            # Get system-wide storage
            total_storage = 1099511627776  # 1TB default
            used_storage = self.db.query(func.sum(File.file_size)).scalar() or 0

        available_storage = total_storage - used_storage
        usage_percentage = (used_storage / total_storage * 100) if total_storage > 0 else 0.0

        # Storage by project. SQLAlchemy 2.x can't auto-figure the join
        # path because Project<->File goes through Sample (no direct FK)
        # — declare each step explicitly.
        project_storage = (
            self.db.query(Project.name, func.sum(File.file_size).label("size"))
            .select_from(Project)
            .join(Sample, Sample.project_id == Project.id)
            .join(File, File.sample_id == Sample.id)
            .group_by(Project.id, Project.name)
            .order_by(desc("size"))
            .limit(10)
            .all()
        )

        by_project = [{"project": row.name, "size": row.size or 0} for row in project_storage]

        # Storage by file type
        file_type_storage = (
            self.db.query(File.file_type, func.sum(File.file_size).label("size")).group_by(File.file_type).all()
        )

        by_file_type = {row.file_type: row.size or 0 for row in file_type_storage}

        # Calculate growth rate (bytes per day) - simplified
        # In production, this would analyze historical data
        growth_rate = None
        projected_full_date = None

        return StorageAnalytics(
            total_storage=total_storage,
            used_storage=used_storage,
            available_storage=available_storage,
            usage_percentage=usage_percentage,
            by_project=by_project,
            by_file_type=by_file_type,
            growth_rate=growth_rate,
            projected_full_date=projected_full_date,
        )

    # ========================================================================
    # Comparative Analytics
    # ========================================================================

    def get_comparative_analysis(
        self, entity_type: str, metric: str, limit: int = 10, user_id: Optional[str] = None
    ) -> ComparativeAnalysis:
        """Get comparative analysis for entities"""

        entities = []

        if entity_type == "project" and metric == "success_rate":
            # Get project success rates
            query = (
                self.db.query(
                    Project.id,
                    Project.name,
                    func.count(PipelineTask.id).label("total_tasks"),
                    func.sum(case((PipelineTask.status == "completed", 1), else_=0)).label("completed"),
                )
                .join(Sample)
                .join(PipelineTask)
            )

            if user_id:
                query = query.filter(Project.owner_id == user_id)

            results = (
                query.group_by(Project.id, Project.name)
                .having(func.count(PipelineTask.id) > 0)
                .order_by(desc("completed"))
                .limit(limit)
                .all()
            )

            for row in results:
                success_rate = (row.completed / row.total_tasks * 100) if row.total_tasks > 0 else 0.0
                entities.append(
                    ComparisonMetric(
                        entity_id=str(row.id), entity_name=row.name, metric_name=metric, value=success_rate, unit="%"
                    )
                )

        elif entity_type == "pipeline" and metric == "avg_duration":
            # Get pipeline average durations
            results = self.get_pipeline_performance(user_id=user_id)
            for perf in results[:limit]:
                if perf.avg_duration:
                    entities.append(
                        ComparisonMetric(
                            entity_id=perf.pipeline_name,
                            entity_name=perf.pipeline_name,
                            metric_name=metric,
                            value=perf.avg_duration,
                            unit="seconds",
                        )
                    )

        # Assign ranks
        entities.sort(key=lambda x: x.value, reverse=True)
        for i, entity in enumerate(entities, 1):
            entity.rank = i

        best = entities[0] if entities else None
        worst = entities[-1] if entities else None
        avg = mean([e.value for e in entities]) if entities else 0.0

        return ComparativeAnalysis(
            entity_type=entity_type,
            metric=metric,
            entities=entities,
            best_performer=best,
            worst_performer=worst,
            average=avg,
            generated_at=datetime.utcnow(),
        )

    # ========================================================================
    # Dashboard Analytics
    # ========================================================================

    def get_dashboard_analytics(self, user_id: Optional[str] = None) -> DashboardAnalytics:
        """Get dashboard analytics summary"""

        # Key metrics
        key_metrics = []

        # Projects metric
        projects_count = self.db.query(func.count(Project.id))
        if user_id:
            projects_count = projects_count.filter(Project.owner_id == user_id)
        projects_total = projects_count.scalar()

        key_metrics.append(
            DashboardMetric(name="Total Projects", value=float(projects_total), unit="projects", trend="stable")
        )

        # Tasks metric
        tasks_query = self.db.query(
            func.count(PipelineTask.id), func.sum(case((PipelineTask.status == "completed", 1), else_=0))
        )
        if user_id:
            tasks_query = tasks_query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        tasks_result = tasks_query.first()
        tasks_total = tasks_result[0] or 0
        tasks_completed = tasks_result[1] or 0
        success_rate = (tasks_completed / tasks_total * 100) if tasks_total > 0 else 0.0

        key_metrics.append(
            DashboardMetric(
                name="Task Success Rate", value=success_rate, unit="%", trend="up" if success_rate > 90 else "down"
            )
        )

        # Storage metric
        storage = self.get_storage_analytics(user_id)
        if storage:
            key_metrics.append(
                DashboardMetric(
                    name="Storage Usage",
                    value=storage.usage_percentage,
                    unit="%",
                    trend="up" if storage.usage_percentage > 80 else "stable",
                )
            )

        # Recent activity (last 10 tasks)
        activity_query = self.db.query(PipelineTask).order_by(desc(PipelineTask.created_at)).limit(10)
        if user_id:
            activity_query = activity_query.join(Sample).join(Project).filter(Project.owner_id == user_id)

        recent_tasks = activity_query.all()
        recent_activity = [
            {
                "type": "task",
                "action": f"{task.name} - {task.status}",
                "timestamp": task.created_at.isoformat(),
                "status": task.status,
            }
            for task in recent_tasks
        ]

        # Alerts
        alerts = []
        if storage and storage.usage_percentage > 80:
            alerts.append(
                {
                    "type": "warning",
                    "message": f"Storage usage at {storage.usage_percentage:.1f}%",
                    "action": "Consider cleaning up old files",
                }
            )

        # Recommendations
        recommendations = []
        if success_rate < 90:
            recommendations.append("Review failed tasks to identify common issues")
        if storage and storage.usage_percentage > 80:
            recommendations.append("Archive or delete unused data to free up space")

        return DashboardAnalytics(
            key_metrics=key_metrics,
            recent_activity=recent_activity,
            alerts=alerts,
            recommendations=recommendations,
            generated_at=datetime.utcnow(),
        )

    # ========================================================================
    # Trend Analysis
    # ========================================================================

    def get_trend_analysis(self, metric: str, time_range: TimeRange, user_id: Optional[str] = None) -> TrendAnalysis:
        """Get trend analysis for a metric"""

        # Get time series data
        time_series = self.get_time_series_data(metric, time_range, user_id=user_id)

        # Simple trend calculation
        if len(time_series.data_points) < 2:
            direction = "stable"
            slope = 0.0
            confidence = 0.0
        else:
            values = [dp.value for dp in time_series.data_points]
            first_half = mean(values[: len(values) // 2])
            second_half = mean(values[len(values) // 2 :])

            if second_half > first_half * 1.1:
                direction = "increasing"
            elif second_half < first_half * 0.9:
                direction = "decreasing"
            else:
                direction = "stable"

            slope = (second_half - first_half) / len(values)
            confidence = 0.8  # Placeholder

        return TrendAnalysis(
            metric=metric,
            time_range=time_range,
            direction=direction,
            slope=slope,
            confidence=confidence,
            data_points=time_series.data_points,
            forecast=None,  # Could implement forecasting here
            anomalies=None,
        )
