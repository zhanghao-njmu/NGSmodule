"""
Analytics schemas for data analysis and visualization
"""
from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any
from datetime import datetime
from enum import Enum


class AnalysisType(str, Enum):
    """Analysis type enumeration"""
    TIME_SERIES = "time_series"
    COMPARISON = "comparison"
    DISTRIBUTION = "distribution"
    CORRELATION = "correlation"
    TREND = "trend"


class TimeRange(str, Enum):
    """Time range enumeration"""
    HOUR = "hour"
    DAY = "day"
    WEEK = "week"
    MONTH = "month"
    YEAR = "year"
    CUSTOM = "custom"


class MetricType(str, Enum):
    """Metric type enumeration"""
    TASKS = "tasks"
    SAMPLES = "samples"
    STORAGE = "storage"
    PERFORMANCE = "performance"
    QUALITY = "quality"


# ============================================================================
# Time Series Analytics
# ============================================================================

class TimeSeriesDataPoint(BaseModel):
    """Single data point in time series"""
    timestamp: datetime
    value: float
    label: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class TimeSeriesData(BaseModel):
    """Time series data"""
    metric: str
    time_range: TimeRange
    data_points: List[TimeSeriesDataPoint]
    total_points: int
    start_time: datetime
    end_time: datetime
    statistics: Optional[Dict[str, float]] = None  # min, max, avg, std


# ============================================================================
# Project Analytics
# ============================================================================

class ProjectPerformanceMetrics(BaseModel):
    """Project performance metrics"""
    project_id: str
    project_name: str
    total_samples: int
    total_tasks: int
    completed_tasks: int
    failed_tasks: int
    success_rate: float
    avg_task_duration: Optional[float] = None  # in seconds
    total_processing_time: Optional[float] = None  # in seconds
    storage_used: int  # in bytes
    created_at: datetime
    last_activity: Optional[datetime] = None


class ProjectComparison(BaseModel):
    """Compare multiple projects"""
    projects: List[ProjectPerformanceMetrics]
    comparison_metric: str
    generated_at: datetime


# ============================================================================
# Sample Analytics
# ============================================================================

class SampleQualityMetrics(BaseModel):
    """Sample quality metrics"""
    sample_id: str
    sample_name: str
    project_id: str
    quality_score: Optional[float] = None  # 0-100
    read_count: Optional[int] = None
    gc_content: Optional[float] = None
    duplication_rate: Optional[float] = None
    q20_percentage: Optional[float] = None
    q30_percentage: Optional[float] = None
    contamination_rate: Optional[float] = None
    status: str
    created_at: datetime


class SampleQualityDistribution(BaseModel):
    """Distribution of sample quality metrics"""
    metric_name: str
    bins: List[Dict[str, Any]]  # [{"range": "0-10", "count": 5}, ...]
    total_samples: int
    statistics: Dict[str, float]  # min, max, avg, median, std


# ============================================================================
# Task Analytics
# ============================================================================

class TaskExecutionMetrics(BaseModel):
    """Task execution metrics"""
    task_id: str
    task_name: str
    pipeline: str
    status: str
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    duration: Optional[float] = None  # in seconds
    cpu_usage: Optional[float] = None
    memory_usage: Optional[float] = None
    exit_code: Optional[int] = None
    retries: int = 0


class PipelinePerformanceMetrics(BaseModel):
    """Pipeline performance metrics"""
    pipeline_name: str
    total_runs: int
    successful_runs: int
    failed_runs: int
    success_rate: float
    avg_duration: Optional[float] = None  # in seconds
    min_duration: Optional[float] = None
    max_duration: Optional[float] = None
    avg_cpu_usage: Optional[float] = None
    avg_memory_usage: Optional[float] = None


class TaskExecutionTrend(BaseModel):
    """Task execution trend over time"""
    period: TimeRange
    data: List[Dict[str, Any]]  # [{"date": "2024-01-01", "completed": 10, "failed": 2}, ...]
    total_periods: int
    generated_at: datetime


# ============================================================================
# Resource Analytics
# ============================================================================

class ResourceUtilization(BaseModel):
    """Resource utilization metrics"""
    timestamp: datetime
    cpu_usage: float  # percentage
    memory_usage: float  # percentage
    disk_usage: float  # percentage
    network_in: Optional[float] = None  # bytes/sec
    network_out: Optional[float] = None  # bytes/sec
    active_tasks: int
    queue_length: int


class ResourceUtilizationTrend(BaseModel):
    """Resource utilization trend"""
    resource_type: str  # cpu, memory, disk
    time_range: TimeRange
    data_points: List[ResourceUtilization]
    peak_usage: float
    avg_usage: float
    generated_at: datetime


class StorageAnalytics(BaseModel):
    """Storage analytics"""
    total_storage: int  # bytes
    used_storage: int  # bytes
    available_storage: int  # bytes
    usage_percentage: float
    by_project: List[Dict[str, Any]]  # [{"project": "P1", "size": 1000}, ...]
    by_file_type: Dict[str, int]  # {"fastq": 1000, "bam": 500}
    growth_rate: Optional[float] = None  # bytes per day
    projected_full_date: Optional[datetime] = None


# ============================================================================
# Correlation Analytics
# ============================================================================

class CorrelationData(BaseModel):
    """Correlation between two metrics"""
    metric_x: str
    metric_y: str
    correlation_coefficient: float  # -1 to 1
    p_value: Optional[float] = None
    data_points: List[Dict[str, float]]  # [{"x": 10, "y": 20}, ...]
    sample_size: int


class MultiMetricCorrelation(BaseModel):
    """Correlation matrix for multiple metrics"""
    metrics: List[str]
    correlation_matrix: List[List[float]]
    generated_at: datetime


# ============================================================================
# Distribution Analytics
# ============================================================================

class DistributionData(BaseModel):
    """Distribution of a metric"""
    metric: str
    bins: List[Dict[str, Any]]  # [{"range": "0-10", "count": 5, "percentage": 25}, ...]
    total_count: int
    statistics: Dict[str, float]  # min, max, mean, median, mode, std, skewness, kurtosis


# ============================================================================
# Comparative Analytics
# ============================================================================

class ComparisonMetric(BaseModel):
    """Single metric for comparison"""
    entity_id: str
    entity_name: str
    metric_name: str
    value: float
    unit: Optional[str] = None
    rank: Optional[int] = None


class ComparativeAnalysis(BaseModel):
    """Comparative analysis results"""
    entity_type: str  # project, sample, pipeline, user
    metric: str
    entities: List[ComparisonMetric]
    best_performer: ComparisonMetric
    worst_performer: ComparisonMetric
    average: float
    generated_at: datetime


# ============================================================================
# Custom Reports
# ============================================================================

class ReportType(str, Enum):
    """Report type enumeration"""
    PROJECT_SUMMARY = "project_summary"
    SAMPLE_QUALITY = "sample_quality"
    PIPELINE_PERFORMANCE = "pipeline_performance"
    RESOURCE_USAGE = "resource_usage"
    CUSTOM = "custom"


class ReportFilter(BaseModel):
    """Report filter parameters"""
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    project_ids: Optional[List[str]] = None
    sample_ids: Optional[List[str]] = None
    pipeline_names: Optional[List[str]] = None
    status: Optional[List[str]] = None


class ReportSection(BaseModel):
    """Single section in a report"""
    title: str
    type: str  # chart, table, text, metrics
    data: Any
    description: Optional[str] = None


class CustomReport(BaseModel):
    """Custom analytics report"""
    report_id: str
    report_type: ReportType
    title: str
    description: Optional[str] = None
    filters: ReportFilter
    sections: List[ReportSection]
    generated_at: datetime
    generated_by: str  # user_id


class ReportCreate(BaseModel):
    """Create custom report request"""
    report_type: ReportType
    title: str
    description: Optional[str] = None
    filters: ReportFilter
    sections: List[Dict[str, Any]]  # Section configuration


# ============================================================================
# Dashboard Analytics
# ============================================================================

class DashboardMetric(BaseModel):
    """Single metric for dashboard"""
    name: str
    value: float
    unit: Optional[str] = None
    change: Optional[float] = None  # percentage change from previous period
    trend: Optional[str] = None  # up, down, stable
    sparkline: Optional[List[float]] = None  # mini chart data


class DashboardAnalytics(BaseModel):
    """Dashboard analytics summary"""
    key_metrics: List[DashboardMetric]
    recent_activity: List[Dict[str, Any]]
    alerts: List[Dict[str, Any]]
    recommendations: List[str]
    generated_at: datetime


# ============================================================================
# Trend Analytics
# ============================================================================

class TrendAnalysis(BaseModel):
    """Trend analysis for a metric"""
    metric: str
    time_range: TimeRange
    direction: str  # increasing, decreasing, stable
    slope: float
    confidence: float  # 0-1
    data_points: List[TimeSeriesDataPoint]
    forecast: Optional[List[TimeSeriesDataPoint]] = None
    anomalies: Optional[List[Dict[str, Any]]] = None


# ============================================================================
# Export/Download
# ============================================================================

class ExportFormat(str, Enum):
    """Export format enumeration"""
    JSON = "json"
    CSV = "csv"
    EXCEL = "excel"
    PDF = "pdf"


class ExportRequest(BaseModel):
    """Export analytics data request"""
    data_type: str  # time_series, comparison, distribution, etc.
    format: ExportFormat
    filters: Optional[ReportFilter] = None
    include_charts: bool = False
