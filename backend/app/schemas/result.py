"""
Result schemas for API request/response
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from uuid import UUID
from datetime import datetime


class ResultBase(BaseModel):
    """Base result schema"""

    result_type: str = Field(..., description="Type of result (qc_report, alignment, quantification, de_analysis)")
    # Read from ORM's `result_metadata` attribute, expose as `metadata` in JSON
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Result metadata", alias="result_metadata")


class ResultResponse(ResultBase):
    """Result response schema"""

    id: UUID
    task_id: UUID
    result_path: Optional[str] = None
    created_at: datetime

    model_config = {"from_attributes": True, "populate_by_name": True}


class ResultListResponse(BaseModel):
    """Result list response schema"""

    results: List[ResultResponse]
    total: int
    skip: int
    limit: int


# Visualization schemas
class QCMetrics(BaseModel):
    """Quality control metrics"""

    total_reads: int
    quality_score: float
    gc_content: float
    duplication_rate: float


class AlignmentStats(BaseModel):
    """Alignment statistics"""

    mapped_reads: int
    mapping_rate: float
    properly_paired: int
    average_coverage: float


class ChartData(BaseModel):
    """Generic chart data structure"""

    type: str = Field(..., description="Chart type (line, bar, pie, scatter, histogram, area)")
    x: Optional[List[Any]] = Field(None, description="X-axis data")
    y: Optional[List[Any]] = Field(None, description="Y-axis data")
    categories: Optional[List[str]] = Field(None, description="Categories for categorical charts")
    values: Optional[List[float]] = Field(None, description="Values for categorical charts")
    labels: Optional[List[str]] = Field(None, description="Labels for data points")
    bins: Optional[int] = Field(None, description="Number of bins for histograms")
    data: Optional[List[float]] = Field(None, description="Raw data for histograms")


class ResultVisualizationData(BaseModel):
    """Visualization data for results"""

    type: str = Field(..., description="Result type")
    metrics: Dict[str, Any] = Field(default_factory=dict, description="Key metrics")
    charts: Dict[str, Any] = Field(default_factory=dict, description="Chart data")
    status: Optional[str] = Field(None, description="Overall status (pass, warning, fail)")
    message: Optional[str] = Field(None, description="Optional message")
    raw_metadata: Optional[Dict[str, Any]] = Field(None, description="Raw metadata if visualization not implemented")
    top_genes: Optional[List[Dict[str, Any]]] = Field(None, description="Top genes for quantification/DE analysis")
    significant_genes: Optional[List[Dict[str, Any]]] = Field(None, description="Significant genes for DE analysis")
