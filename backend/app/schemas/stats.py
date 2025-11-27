"""
Statistics schemas for API responses
"""
from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime


class ProjectStats(BaseModel):
    """Project statistics"""
    total: int = Field(..., description="Total number of projects")
    active: int = Field(..., description="Number of active projects")
    completed: int = Field(..., description="Number of completed projects")
    failed: int = Field(..., description="Number of failed projects")

    class Config:
        from_attributes = True


class SampleStats(BaseModel):
    """Sample statistics"""
    total: int = Field(..., description="Total number of samples")
    processing: int = Field(..., description="Number of samples being processed")
    completed: int = Field(..., description="Number of completed samples")
    failed: int = Field(..., description="Number of failed samples")

    class Config:
        from_attributes = True


class TaskStats(BaseModel):
    """Task statistics"""
    total: int = Field(..., description="Total number of tasks")
    pending: int = Field(..., description="Number of pending tasks")
    running: int = Field(..., description="Number of running tasks")
    completed: int = Field(..., description="Number of completed tasks")
    failed: int = Field(..., description="Number of failed tasks")
    cancelled: int = Field(..., description="Number of cancelled tasks")

    class Config:
        from_attributes = True


class FileStats(BaseModel):
    """File statistics"""
    total: int = Field(..., description="Total number of files")
    total_size: int = Field(..., description="Total size in bytes")
    by_type: Dict[str, int] = Field(default_factory=dict, description="Files grouped by type")

    class Config:
        from_attributes = True


class StorageStats(BaseModel):
    """Storage usage statistics"""
    total_space: int = Field(..., description="Total storage space in bytes")
    used_space: int = Field(..., description="Used storage space in bytes")
    available_space: int = Field(..., description="Available storage space in bytes")
    usage_percentage: float = Field(..., description="Storage usage percentage")
    files_count: int = Field(..., description="Total number of files")
    projects_space: int = Field(..., description="Space used by projects")
    samples_space: int = Field(..., description="Space used by samples")
    results_space: int = Field(..., description="Space used by results")

    class Config:
        from_attributes = True


class UserActivityStats(BaseModel):
    """User activity statistics"""
    total_users: int = Field(..., description="Total number of users")
    active_users_today: int = Field(..., description="Active users today")
    active_users_week: int = Field(..., description="Active users this week")
    active_users_month: int = Field(..., description="Active users this month")
    new_users_week: int = Field(..., description="New users this week")
    new_users_month: int = Field(..., description="New users this month")

    class Config:
        from_attributes = True


class PipelineStats(BaseModel):
    """Pipeline execution statistics"""
    total_runs: int = Field(..., description="Total pipeline runs")
    success_rate: float = Field(..., description="Success rate percentage")
    average_duration: float = Field(..., description="Average duration in seconds")
    most_used: List[Dict[str, Any]] = Field(default_factory=list, description="Most used pipelines")

    class Config:
        from_attributes = True


class SystemStats(BaseModel):
    """Overall system statistics"""
    uptime: float = Field(..., description="System uptime in seconds")
    api_calls_today: int = Field(..., description="API calls today")
    api_calls_hour: int = Field(..., description="API calls in the last hour")
    active_connections: int = Field(..., description="Active WebSocket connections")
    queue_size: int = Field(..., description="Task queue size")

    class Config:
        from_attributes = True


class StatsSummary(BaseModel):
    """Comprehensive statistics summary for dashboard"""
    projects: ProjectStats
    samples: SampleStats
    tasks: TaskStats
    files: FileStats
    storage: StorageStats
    users: Optional[UserActivityStats] = None
    pipelines: Optional[PipelineStats] = None
    system: Optional[SystemStats] = None
    generated_at: datetime = Field(default_factory=datetime.utcnow)

    class Config:
        from_attributes = True


class TrendData(BaseModel):
    """Trend data for time-series analysis"""
    metric: str = Field(..., description="Metric name")
    period: str = Field(..., description="Time period (day, week, month)")
    data_points: List[Dict[str, Any]] = Field(..., description="Time-series data points")

    class Config:
        from_attributes = True


class QuickStats(BaseModel):
    """Quick statistics for header/widget display"""
    total_projects: int
    total_samples: int
    running_tasks: int
    storage_used: int
    storage_quota: int

    class Config:
        from_attributes = True
