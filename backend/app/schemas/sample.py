"""
Sample schemas for API request/response
"""
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any, List
from datetime import datetime
from uuid import UUID


class SampleBase(BaseModel):
    """Base sample schema"""
    sample_id: str = Field(..., min_length=1, max_length=100)
    run_id: Optional[str] = Field(None, max_length=100)
    group_name: Optional[str] = Field(None, max_length=50)
    layout: Optional[str] = Field(None, description="Layout: PE (paired-end) or SE (single-end)")
    batch_id: Optional[str] = Field(None, max_length=50)


class SampleCreate(SampleBase):
    """Schema for creating a new sample"""
    project_id: UUID
    metadata: Optional[Dict[str, Any]] = Field(default_factory=dict)


class SampleBatchCreate(BaseModel):
    """Schema for batch creating samples"""
    project_id: UUID
    samples: List[SampleBase]


class SampleUpdate(BaseModel):
    """Schema for updating sample"""
    sample_id: Optional[str] = Field(None, min_length=1, max_length=100)
    run_id: Optional[str] = None
    group_name: Optional[str] = None
    layout: Optional[str] = None
    batch_id: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class SampleResponse(SampleBase):
    """Schema for sample response"""
    id: UUID
    project_id: UUID
    # Read from ORM's sample_metadata attribute, expose as `metadata` in JSON
    metadata: Dict[str, Any] = Field(default_factory=dict, alias="sample_metadata")
    created_at: datetime

    # Computed fields
    file_count: Optional[int] = None

    class Config:
        from_attributes = True
        populate_by_name = True


class SampleListResponse(BaseModel):
    """Schema for sample list response"""
    total: int
    items: List[SampleResponse]
    page: int
    page_size: int


class SampleImportFromCSV(BaseModel):
    """Schema for importing samples from CSV"""
    project_id: UUID
    csv_content: str = Field(..., description="CSV content with headers: sample_id,run_id,group_name,layout,batch_id")
