"""
Pipeline Template schemas
"""

from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel, Field


class PipelineTemplateBase(BaseModel):
    """Base pipeline template schema"""

    name: str = Field(..., min_length=1, max_length=100)
    display_name: str = Field(..., min_length=1, max_length=200)
    description: Optional[str] = None
    category: str = Field(..., min_length=1, max_length=50)
    script_name: str = Field(..., min_length=1, max_length=200)
    script_path: Optional[str] = None
    default_params: Dict[str, Any] = Field(default_factory=dict)
    param_schema: Dict[str, Any] = Field(default_factory=dict)
    estimated_time: Optional[str] = None
    min_memory_gb: Optional[str] = None
    min_cpu_cores: Optional[str] = None
    sort_order: str = Field(default="0")
    tags: List[str] = Field(default_factory=list)


class PipelineTemplateCreate(PipelineTemplateBase):
    """Schema for creating pipeline template"""


class PipelineTemplateUpdate(BaseModel):
    """Schema for updating pipeline template"""

    display_name: Optional[str] = None
    description: Optional[str] = None
    category: Optional[str] = None
    script_path: Optional[str] = None
    default_params: Optional[Dict[str, Any]] = None
    param_schema: Optional[Dict[str, Any]] = None
    estimated_time: Optional[str] = None
    min_memory_gb: Optional[str] = None
    min_cpu_cores: Optional[str] = None
    is_active: Optional[bool] = None
    sort_order: Optional[str] = None
    tags: Optional[List[str]] = None


class PipelineTemplateResponse(PipelineTemplateBase):
    """Schema for pipeline template response"""

    id: UUID
    is_active: bool
    is_builtin: bool

    model_config = {"from_attributes": True}


class PipelineTemplateListResponse(BaseModel):
    """Schema for pipeline template list response"""

    total: int
    items: List[PipelineTemplateResponse]


class PipelineExecuteRequest(BaseModel):
    """Schema for executing a pipeline"""

    template_id: UUID = Field(..., description="Pipeline template ID")
    task_name: str = Field(..., min_length=1, max_length=100, description="Task name")
    project_id: UUID = Field(..., description="Project ID")
    sample_ids: List[UUID] = Field(default_factory=list, description="Sample IDs to process")
    parameters: Dict[str, Any] = Field(default_factory=dict, description="Pipeline parameters")


class PipelineTemplateCategory(BaseModel):
    """Schema for pipeline categories"""

    category: str
    count: int
    templates: List[str]


class PipelineBatchExecuteRequest(BaseModel):
    """Schema for batch executing a pipeline on multiple samples"""

    template_id: UUID = Field(..., description="Pipeline template ID")
    project_id: UUID = Field(..., description="Project ID")
    sample_ids: List[UUID] = Field(..., min_length=1, description="Sample IDs to process (one task per sample)")
    task_name_prefix: str = Field(
        ..., min_length=1, max_length=80, description="Task name prefix (will append sample name)"
    )
    parameters: Dict[str, Any] = Field(default_factory=dict, description="Pipeline parameters (same for all samples)")


class PipelineBatchExecuteResponse(BaseModel):
    """Schema for batch execution response"""

    total_tasks: int
    created_tasks: List[UUID]
    failed_samples: List[Dict[str, str]] = Field(default_factory=list)


class ParameterRecommendationResponse(BaseModel):
    """Schema for parameter recommendation response"""

    recommended_params: Dict[str, Any]
    confidence_score: float = Field(..., ge=0.0, le=1.0, description="Confidence score (0-1)")
    based_on_tasks: int = Field(..., description="Number of historical tasks analyzed")
    explanation: str = Field(..., description="Explanation of recommendation")
