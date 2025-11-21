"""
File schemas for API request/response
"""
from pydantic import BaseModel, Field
from typing import Optional
from datetime import datetime
from uuid import UUID


class FileBase(BaseModel):
    """Base file schema"""
    filename: str = Field(..., min_length=1, max_length=255)
    file_type: Optional[str] = Field(None, max_length=20)


class FileCreate(BaseModel):
    """Schema for creating file record"""
    sample_id: UUID
    filename: str
    file_path: str
    file_type: Optional[str] = None
    file_size: Optional[int] = None
    md5_checksum: Optional[str] = None


class FileUpdate(BaseModel):
    """Schema for updating file"""
    upload_status: Optional[str] = Field(
        None,
        description="Status: pending, uploading, completed, failed"
    )
    file_size: Optional[int] = None
    md5_checksum: Optional[str] = None


class FileResponse(FileBase):
    """Schema for file response"""
    id: UUID
    sample_id: UUID
    file_path: str
    file_size: Optional[int] = None
    md5_checksum: Optional[str] = None
    upload_status: str
    created_at: datetime

    class Config:
        from_attributes = True


class FileListResponse(BaseModel):
    """Schema for file list response"""
    total: int
    items: list[FileResponse]


class FileUploadInitiate(BaseModel):
    """Schema for initiating file upload"""
    sample_id: UUID
    filename: str
    file_size: int
    file_type: Optional[str] = None
    md5_checksum: Optional[str] = None


class FileUploadComplete(BaseModel):
    """Schema for completing file upload"""
    file_id: UUID
    md5_checksum: str
