"""
Statistics API endpoints
"""
from fastapi import APIRouter, Depends, Query
from sqlalchemy.orm import Session
from typing import Optional

from app.core.deps import get_db, get_current_user, get_current_active_user
from app.core.permissions import require_role
from app.models.user import User
from app.services.stats_service import StatsService
from app.schemas.stats import (
    StatsSummary,
    ProjectStats,
    SampleStats,
    TaskStats,
    FileStats,
    StorageStats,
    UserActivityStats,
    PipelineStats,
    SystemStats,
    QuickStats,
    TrendData
)

router = APIRouter()


@router.get("/summary", response_model=StatsSummary)
async def get_stats_summary(
    include_system: bool = Query(False, description="Include system-level statistics"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get comprehensive statistics summary

    - **include_system**: Include system-level stats (admin only)
    """
    # Only admin can request system stats
    if include_system and current_user.role != 'admin':
        include_system = False

    service = StatsService(db)
    return service.get_summary(include_system=include_system)


@router.get("/projects", response_model=ProjectStats)
async def get_project_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get project statistics"""
    service = StatsService(db)
    return service.get_project_stats()


@router.get("/samples", response_model=SampleStats)
async def get_sample_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get sample statistics"""
    service = StatsService(db)
    return service.get_sample_stats()


@router.get("/tasks", response_model=TaskStats)
async def get_task_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get task statistics"""
    service = StatsService(db)
    return service.get_task_stats()


@router.get("/files", response_model=FileStats)
async def get_file_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get file statistics"""
    service = StatsService(db)
    return service.get_file_stats()


@router.get("/storage", response_model=StorageStats)
async def get_storage_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get storage usage statistics"""
    service = StatsService(db)
    return service.get_storage_stats()


@router.get("/users", response_model=UserActivityStats)
@require_role("admin")
async def get_user_activity_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get user activity statistics (admin only)
    """
    service = StatsService(db)
    return service.get_user_activity_stats()


@router.get("/pipelines", response_model=PipelineStats)
async def get_pipeline_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """Get pipeline execution statistics"""
    service = StatsService(db)
    return service.get_pipeline_stats()


@router.get("/system", response_model=SystemStats)
@require_role("admin")
async def get_system_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get system-level statistics (admin only)
    """
    service = StatsService(db)
    return service.get_system_stats()


@router.get("/quick", response_model=QuickStats)
async def get_quick_stats(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get quick statistics for header/widget display

    Returns essential metrics for quick display
    """
    service = StatsService(db)
    return service.get_quick_stats(user_id=str(current_user.id))


@router.get("/trends/{metric}", response_model=TrendData)
async def get_trend_data(
    metric: str,
    period: str = Query("week", pattern="^(day|week|month)$"),
    days: int = Query(7, ge=1, le=90),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user)
):
    """
    Get trend data for time-series analysis

    - **metric**: Metric name (projects, samples, tasks, etc.)
    - **period**: Time period (day, week, month)
    - **days**: Number of days to include (1-90)
    """
    service = StatsService(db)
    return service.get_trend_data(metric=metric, period=period, days=days)
