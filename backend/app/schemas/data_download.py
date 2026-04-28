"""
Data download schemas (API request/response).
"""
from datetime import datetime
from typing import Optional, Literal
from uuid import UUID

from pydantic import BaseModel, Field


DownloadStatus = Literal["pending", "running", "completed", "failed", "cancelled"]


# ---------------- session ----------------

class SessionLogin(BaseModel):
    """Open a vendor session (e.g. lcbio `run email password`)."""
    vendor: str = Field(..., description="Vendor id, see /data-downloads/vendors")
    email: str
    password: str = Field(..., min_length=1)


class SessionStatus(BaseModel):
    """Current state of a vendor session (i.e. is its long-running daemon up)."""
    vendor: str
    active: bool
    pid: Optional[int] = None


# ---------------- jobs ----------------

class DownloadJobCreate(BaseModel):
    """Start a new download. Vendor session must already be active."""
    vendor: str = Field(..., description="Vendor id, see /data-downloads/vendors")
    source_path: str = Field(..., min_length=1, description="Vendor-side path; for lcbio this is the obs_path copied from web UI")
    dest_path: str = Field(..., min_length=1, description="Local destination directory")


class DownloadJobResponse(BaseModel):
    """Single download job state."""
    id: UUID
    vendor: str
    source_path: str
    dest_path: str
    log_path: Optional[str] = None
    status: DownloadStatus
    progress_pct: float
    bytes_downloaded: Optional[int] = None
    file_size: Optional[int] = None
    error_message: Optional[str] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    created_at: datetime

    class Config:
        from_attributes = True


class DownloadJobListResponse(BaseModel):
    total: int
    items: list[DownloadJobResponse]
