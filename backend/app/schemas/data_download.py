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
    """Open a vendor session (e.g. lcbio `run email password`).

    Either pass `email` + `password`, or `credential_id` referring to a
    saved entry from POST /api/v1/vendor-credentials.
    """
    vendor: str = Field(..., description="Vendor id, see /data-downloads/vendors")
    email: Optional[str] = None
    password: Optional[str] = Field(default=None, min_length=1)
    credential_id: Optional[UUID] = Field(
        default=None,
        description="If set, decrypt and use the saved credential instead of the inline email/password.",
    )


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
    auto_register: bool = Field(
        default=False,
        description="If true and the delivery is a tar archive, the worker will "
                    "extract it on completion and create an NGSmodule Project "
                    "scaffold pointing at the extracted files (samples are NOT "
                    "auto-created; user reviews them via UI).",
    )
    project_name: Optional[str] = Field(
        default=None,
        max_length=100,
        description="Project name when auto_register is true. If omitted, "
                    "derived from the first segment of source_path.",
    )


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
    project_id: Optional[UUID] = None  # set when auto_register triggered

    model_config = {"from_attributes": True}
class DownloadJobListResponse(BaseModel):
    total: int
    items: list[DownloadJobResponse]
