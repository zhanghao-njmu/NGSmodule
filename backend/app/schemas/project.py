"""
Project schemas for API request/response
"""
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any
from datetime import datetime
from uuid import UUID


class ProjectBase(BaseModel):
    """Base project schema"""
    name: str = Field(..., min_length=1, max_length=100)
    description: Optional[str] = None
    project_type: Optional[str] = Field(
        None,
        description="Project type: rna-seq, dna-seq, sc-rna-seq, atac-seq, chip-seq, etc."
    )


class ProjectCreate(ProjectBase):
    """Schema for creating a new project"""
    config: Optional[Dict[str, Any]] = Field(default_factory=dict)


class ProjectUpdate(BaseModel):
    """Schema for updating project"""
    name: Optional[str] = Field(None, min_length=1, max_length=100)
    description: Optional[str] = None
    project_type: Optional[str] = None
    status: Optional[str] = Field(None, description="Status: active, archived, deleted")
    config: Optional[Dict[str, Any]] = None


class ProjectResponse(ProjectBase):
    """Schema for project response"""
    id: UUID
    user_id: UUID
    status: str
    config: Dict[str, Any]
    created_at: datetime
    updated_at: datetime

    # Computed fields
    sample_count: Optional[int] = None
    task_count: Optional[int] = None

    model_config = {"from_attributes": True}
class ProjectListResponse(BaseModel):
    """Schema for project list response"""
    total: int
    items: list[ProjectResponse]
    page: int
    page_size: int


class ProjectStats(BaseModel):
    """Project statistics"""
    total_projects: int
    active_projects: int
    archived_projects: int
    total_samples: int
    total_tasks: int
    running_tasks: int
    completed_tasks: int
    failed_tasks: int
