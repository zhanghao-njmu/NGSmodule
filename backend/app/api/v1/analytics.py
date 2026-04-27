"""
Analytics API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session
from typing import List, Optional
from datetime import datetime

from app.core.database import get_db
from app.core.deps import get_current_user, get_current_admin as get_current_admin_user
from app.models.user import User
from app.services.analytics_service import AnalyticsService
from app.schemas.analytics import *


router = APIRouter()


# ============================================================================
# Time Series Analytics
# ============================================================================

@router.get("/timeseries/{metric}", response_model=TimeSeriesData)
async def get_time_series(
    metric: str,
    time_range: TimeRange = Query(TimeRange.WEEK),
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get time series data for a metric

    **Metrics:**
    - tasks: Task count over time
    - samples: Sample count over time
    - storage: Storage usage over time

    **Time Ranges:**
    - hour, day, week, month, year, custom
    """
    service = AnalyticsService(db)

    # Non-admin users see only their own data
    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_time_series_data(
        metric=metric,
        time_range=time_range,
        start_date=start_date,
        end_date=end_date,
        user_id=user_id
    )


# ============================================================================
# Project Analytics
# ============================================================================

@router.get("/projects/{project_id}/performance", response_model=ProjectPerformanceMetrics)
async def get_project_performance(
    project_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get performance metrics for a specific project

    Returns:
    - Total samples, tasks, completion rates
    - Average task duration
    - Storage usage
    - Success rate
    """
    service = AnalyticsService(db)

    # Non-admin users can only see their own projects
    user_id = None if current_user.is_admin else str(current_user.id)

    metrics = service.get_project_performance(project_id, user_id)
    if not metrics:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found or access denied"
        )

    return metrics


@router.post("/projects/compare", response_model=ProjectComparison)
async def compare_projects(
    project_ids: List[str],
    metric: str = Query("success_rate", description="Metric to compare"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Compare multiple projects based on a metric

    **Metrics:**
    - success_rate: Task success rate
    - total_tasks: Total number of tasks
    - storage_used: Storage consumption
    - avg_task_duration: Average task execution time
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.compare_projects(project_ids, metric, user_id)


# ============================================================================
# Sample Quality Analytics
# ============================================================================

@router.get("/samples/quality/distribution", response_model=SampleQualityDistribution)
async def get_sample_quality_distribution(
    metric_name: str = Query("quality_score", description="Quality metric name"),
    project_id: Optional[str] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get distribution of sample quality metrics

    **Metrics:**
    - quality_score: Overall quality score (0-100)
    - q20_percentage: Percentage of bases with Q20+
    - q30_percentage: Percentage of bases with Q30+
    - gc_content: GC content percentage
    - duplication_rate: Read duplication rate
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_sample_quality_distribution(
        metric_name=metric_name,
        project_id=project_id,
        user_id=user_id
    )


# ============================================================================
# Task Execution Analytics
# ============================================================================

@router.get("/pipelines/performance", response_model=List[PipelinePerformanceMetrics])
async def get_pipeline_performance(
    pipeline_name: Optional[str] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get performance metrics for pipelines

    Returns metrics for all pipelines or a specific pipeline:
    - Total runs, success/failure counts
    - Success rate
    - Average, min, max execution duration
    - Resource usage (if available)
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_pipeline_performance(pipeline_name, user_id)


@router.get("/tasks/execution-trend", response_model=TaskExecutionTrend)
async def get_task_execution_trend(
    time_range: TimeRange = Query(TimeRange.WEEK),
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get task execution trend over time

    Returns daily breakdown of:
    - Total tasks created
    - Completed tasks
    - Failed tasks
    - Running tasks
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_task_execution_trend(
        time_range=time_range,
        start_date=start_date,
        end_date=end_date,
        user_id=user_id
    )


# ============================================================================
# Resource Analytics
# ============================================================================

@router.get("/storage", response_model=StorageAnalytics)
async def get_storage_analytics(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get storage analytics

    Returns:
    - Total, used, available storage
    - Usage percentage
    - Storage by project
    - Storage by file type
    - Growth rate and projections
    """
    service = AnalyticsService(db)

    # Non-admin users see only their own storage
    user_id = None if current_user.is_admin else str(current_user.id)

    analytics = service.get_storage_analytics(user_id)
    if not analytics:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Storage analytics not available"
        )

    return analytics


# ============================================================================
# Comparative Analytics
# ============================================================================

@router.get("/compare/{entity_type}", response_model=ComparativeAnalysis)
async def get_comparative_analysis(
    entity_type: str,
    metric: str = Query(..., description="Metric to compare"),
    limit: int = Query(10, ge=1, le=100, description="Number of entities to compare"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get comparative analysis for entities

    **Entity Types:**
    - project: Compare projects
    - pipeline: Compare pipeline types
    - user: Compare users (admin only)

    **Metrics:**
    - success_rate: Task success rate
    - avg_duration: Average execution time
    - total_runs: Total number of runs
    - storage_used: Storage consumption
    """
    service = AnalyticsService(db)

    # Only admins can compare users
    if entity_type == "user" and not current_user.is_admin:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Only administrators can compare users"
        )

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_comparative_analysis(
        entity_type=entity_type,
        metric=metric,
        limit=limit,
        user_id=user_id
    )


# ============================================================================
# Dashboard Analytics
# ============================================================================

@router.get("/dashboard", response_model=DashboardAnalytics)
async def get_dashboard_analytics(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get dashboard analytics summary

    Returns:
    - Key metrics with trends
    - Recent activity
    - Alerts and warnings
    - Recommendations

    Optimized for dashboard display with quick load times.
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_dashboard_analytics(user_id)


# ============================================================================
# Trend Analysis
# ============================================================================

@router.get("/trends/{metric}", response_model=TrendAnalysis)
async def get_trend_analysis(
    metric: str,
    time_range: TimeRange = Query(TimeRange.MONTH),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get trend analysis for a metric

    Returns:
    - Trend direction (increasing, decreasing, stable)
    - Slope and confidence
    - Data points
    - Forecast (if available)
    - Anomalies detected

    **Metrics:**
    - tasks: Task creation trend
    - samples: Sample processing trend
    - storage: Storage growth trend
    - success_rate: Success rate trend
    """
    service = AnalyticsService(db)

    user_id = None if current_user.is_admin else str(current_user.id)

    return service.get_trend_analysis(
        metric=metric,
        time_range=time_range,
        user_id=user_id
    )


# ============================================================================
# Export Analytics
# ============================================================================

@router.post("/export")
async def export_analytics(
    export_request: ExportRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Export analytics data in various formats

    **Formats:**
    - json: JSON format
    - csv: CSV format
    - excel: Excel spreadsheet
    - pdf: PDF report (with charts if requested)

    Returns a download link or file stream.

    Note: This is a placeholder endpoint. Full implementation would
    include actual file generation and download mechanisms.
    """
    # Placeholder response
    return {
        "status": "success",
        "message": f"Export request received for {export_request.data_type} in {export_request.format} format",
        "download_url": "/downloads/analytics_export_123.json",
        "expires_at": datetime.utcnow().isoformat()
    }


# ============================================================================
# Health Check
# ============================================================================

@router.get("/health")
async def analytics_health_check():
    """
    Health check endpoint for analytics service
    """
    return {
        "status": "healthy",
        "service": "analytics",
        "timestamp": datetime.utcnow().isoformat()
    }
