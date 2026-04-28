"""
Data download API endpoints.

Bridges the frontend to vendor-specific download tools (联川 / 诺禾致源 / ...)
through the unified `DataDownloadService`.
"""

from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.schemas.common import MessageResponse
from app.schemas.data_download import (
    DownloadJobCreate,
    DownloadJobListResponse,
    DownloadJobResponse,
    SessionLogin,
    SessionStatus,
)
from app.services.data_download_service import DataDownloadService
from app.services.data_provider.factory import list_vendors

router = APIRouter()


def get_service(db: Session = Depends(get_db)) -> DataDownloadService:
    return DataDownloadService(db)


# ---------------- vendor list ----------------


@router.get("/vendors", response_model=list[str])
async def supported_vendors():
    """List vendor IDs the backend has adapters for."""
    return list_vendors()


# ---------------- session ----------------


@router.get("/sessions/{vendor}", response_model=SessionStatus)
async def get_session(
    vendor: str,
    _: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Whether the vendor's daemon is running. Process-level shared."""
    info = service.session_status(vendor)
    return SessionStatus(vendor=vendor, active=info.active, pid=info.pid)


@router.post("/sessions", response_model=SessionStatus, status_code=status.HTTP_201_CREATED)
async def open_session(
    payload: SessionLogin,
    user: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Start the vendor daemon and authenticate.

    Either inline email+password or `credential_id` from a saved entry.
    """
    if payload.credential_id is not None:
        info = service.login_with_credential(payload.vendor, payload.credential_id, user.id)
    elif payload.email and payload.password:
        info = service.login(payload.vendor, payload.email, payload.password)
    else:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Provide either credential_id or both email and password",
        )
    return SessionStatus(vendor=payload.vendor, active=info.active, pid=info.pid)


@router.delete("/sessions/{vendor}", response_model=MessageResponse)
async def close_session(
    vendor: str,
    _: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Kill the vendor daemon (also stops any in-flight downloads)."""
    service.logout(vendor)
    return MessageResponse(message=f"{vendor} session closed")


# ---------------- jobs ----------------


@router.post("/jobs", response_model=DownloadJobResponse, status_code=status.HTTP_201_CREATED)
async def create_job(
    payload: DownloadJobCreate,
    user: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Trigger a vendor download. Vendor session must be active."""
    job = service.create_job(
        user_id=user.id,
        vendor=payload.vendor,
        source_path=payload.source_path,
        dest_path=payload.dest_path,
        auto_register=payload.auto_register,
        project_name=payload.project_name,
    )
    return job


@router.get("/jobs", response_model=DownloadJobListResponse)
async def list_jobs(
    user: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """List the current user's download jobs (newest first)."""
    jobs = service.list_jobs(user.id)
    # Refresh progress for in-flight jobs so polling clients see fresh data.
    refreshed = [service.refresh_progress(j) if j.status == "running" else j for j in jobs]
    return DownloadJobListResponse(total=len(refreshed), items=refreshed)


@router.get("/jobs/{job_id}", response_model=DownloadJobResponse)
async def get_job(
    job_id: UUID,
    user: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Get a single job's latest state (re-parses the vendor log)."""
    job = service.get_by_id_or_raise(job_id, user.id)
    return service.refresh_progress(job)


@router.delete("/jobs/{job_id}", response_model=DownloadJobResponse)
async def cancel_job(
    job_id: UUID,
    user: User = Depends(get_current_user),
    service: DataDownloadService = Depends(get_service),
):
    """Cancel a running job (lcbio: kills the whole daemon)."""
    job = service.get_by_id_or_raise(job_id, user.id)
    return service.cancel(job)
