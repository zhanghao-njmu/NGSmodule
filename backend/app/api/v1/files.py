"""
File management API endpoints (Refactored with Service Layer)

This is the refactored version using FileService for business logic.
Routers are now thin HTTP adapters that delegate to the service layer.
"""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends
from fastapi import File as FileUpload
from fastapi import Query, Request, UploadFile, status
from fastapi.responses import RedirectResponse
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_user
from app.core.rate_limit import RateLimits, limiter
from app.models.user import User
from app.schemas.file import (
    FileListResponse,
    FileResponse,
)
from app.services.file_service import FileService

router = APIRouter()


# ============= DEPENDENCY INJECTION =============


def get_file_service(db: Session = Depends(get_db)) -> FileService:
    """Dependency to get FileService instance"""
    return FileService(db)


# ============= ENDPOINTS =============


@router.get("", response_model=FileListResponse)
async def list_files(
    sample_id: Optional[UUID] = Query(None, description="Filter by sample"),
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    file_type: Optional[str] = Query(None, description="Filter by file type"),
    current_user: User = Depends(get_current_user),
    service: FileService = Depends(get_file_service),
):
    """
    List files with optional filters
    """
    files = service.list_files(user_id=current_user.id, sample_id=sample_id, project_id=project_id, file_type=file_type)

    return FileListResponse(total=len(files), items=files)


@router.post("/upload", response_model=FileResponse, status_code=status.HTTP_201_CREATED)
@limiter.limit(RateLimits.FILE_UPLOAD)
async def upload_file(
    request: Request,
    sample_id: UUID,
    file: UploadFile = FileUpload(...),
    current_user: User = Depends(get_current_user),
    service: FileService = Depends(get_file_service),
):
    """
    Upload a file and associate with sample

    - **sample_id**: Sample ID to associate file with
    - **file**: File to upload (supports .fastq, .fastq.gz, .bam, .sam, .vcf, etc.)

    Maximum file size: 50GB
    """
    return await service.upload(user_id=current_user.id, sample_id=sample_id, file=file)


@router.get("/{file_id}", response_model=FileResponse)
async def get_file(
    file_id: UUID, current_user: User = Depends(get_current_user), service: FileService = Depends(get_file_service)
):
    """
    Get file information by ID
    """
    return service.get_by_id_or_raise(file_id, current_user.id)


@router.get("/{file_id}/download")
@limiter.limit(RateLimits.FILE_DOWNLOAD)
async def download_file(
    request: Request,
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    service: FileService = Depends(get_file_service),
):
    """
    Download file

    Returns a redirect to presigned URL with the file content
    """
    presigned_url = service.get_download_url(file_id, current_user.id)
    return RedirectResponse(url=presigned_url)


@router.delete("/{file_id}", status_code=status.HTTP_204_NO_CONTENT)
@limiter.limit(RateLimits.FILE_DELETE)
async def delete_file(
    request: Request,
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    service: FileService = Depends(get_file_service),
):
    """
    Delete file

    Removes file from storage and updates user storage quota
    """
    await service.delete(file_id, current_user.id)
    return None
