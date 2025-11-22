"""
File management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query, UploadFile, File as FileUpload, Request
from fastapi.responses import StreamingResponse
from sqlalchemy.orm import Session
from typing import Optional
from uuid import UUID
import os

from app.core.database import get_db
from app.core.deps import get_current_user
from app.core.rate_limit import limiter, RateLimits
from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.file import File
from app.schemas.file import (
    FileCreate,
    FileUpdate,
    FileResponse,
    FileListResponse,
    FileUploadInitiate,
    FileUploadComplete,
)
from app.schemas.common import MessageResponse
from app.services.storage import storage_service
from app.core.config import settings
import hashlib

router = APIRouter()


@router.get("", response_model=FileListResponse)
async def list_files(
    sample_id: Optional[UUID] = Query(None, description="Filter by sample"),
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    file_type: Optional[str] = Query(None, description="Filter by file type"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    List files with optional filters
    """
    query = db.query(File).join(Sample).join(Project).filter(
        Project.user_id == current_user.id
    )

    # Apply filters
    if sample_id:
        query = query.filter(File.sample_id == sample_id)

    if project_id:
        query = query.filter(Project.id == project_id)

    if file_type:
        query = query.filter(File.file_type == file_type)

    files = query.order_by(File.created_at.desc()).all()

    return FileListResponse(
        total=len(files),
        items=files
    )


@router.post("/upload", response_model=FileResponse, status_code=status.HTTP_201_CREATED)
@limiter.limit(RateLimits.FILE_UPLOAD)
async def upload_file(
    request: Request,
    sample_id: UUID,
    file: UploadFile = FileUpload(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Upload a file and associate with sample

    - **sample_id**: Sample ID to associate file with
    - **file**: File to upload (supports .fastq, .fastq.gz, .bam, .sam, .vcf, etc.)

    Maximum file size: 50GB
    """
    # Verify sample belongs to user
    sample = db.query(Sample).join(Project).filter(
        Sample.id == sample_id,
        Project.user_id == current_user.id
    ).first()

    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Sample not found"
        )

    # Check file extension
    file_ext = os.path.splitext(file.filename)[1].lower()
    allowed_extensions = settings.ALLOWED_EXTENSIONS

    if file_ext not in allowed_extensions:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"File type '{file_ext}' not allowed. Allowed: {', '.join(allowed_extensions)}"
        )

    # Check storage quota
    file_size = 0
    try:
        file_contents = await file.read()
        file_size = len(file_contents)

        if current_user.storage_used + file_size > current_user.storage_quota:
            raise HTTPException(
                status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE,
                detail="Storage quota exceeded"
            )

        # Reset file pointer
        await file.seek(0)

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Error reading file: {str(e)}"
        )

    # Save to storage
    try:
        import io
        file_stream = io.BytesIO(file_contents)
        object_path = await storage_service.save_upload_file(
            file=file_stream,
            user_id=str(current_user.id),
            filename=file.filename,
            project_id=str(sample.project_id)
        )

        # Calculate MD5
        md5_hash = hashlib.md5(file_contents).hexdigest()

        # Create file record
        file_record = File(
            sample_id=sample_id,
            filename=file.filename,
            file_path=object_path,
            file_type=file_ext.lstrip('.'),
            file_size=file_size,
            md5_checksum=md5_hash,
            upload_status="completed"
        )

        db.add(file_record)

        # Update user storage
        current_user.storage_used += file_size

        db.commit()
        db.refresh(file_record)

        return file_record

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error uploading file: {str(e)}"
        )


@router.get("/{file_id}", response_model=FileResponse)
async def get_file(
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get file information by ID
    """
    file = db.query(File).join(Sample).join(Project).filter(
        File.id == file_id,
        Project.user_id == current_user.id
    ).first()

    if not file:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    return file


@router.get("/{file_id}/download")
@limiter.limit(RateLimits.FILE_DOWNLOAD)
async def download_file(
    request: Request,
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Download file

    Returns a streaming response with the file content
    """
    file = db.query(File).join(Sample).join(Project).filter(
        File.id == file_id,
        Project.user_id == current_user.id
    ).first()

    if not file:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    try:
        # Get presigned URL (redirect to MinIO)
        presigned_url = storage_service.get_presigned_url(file.file_path)

        # Return redirect to presigned URL
        from fastapi.responses import RedirectResponse
        return RedirectResponse(url=presigned_url)

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error downloading file: {str(e)}"
        )


@router.delete("/{file_id}", status_code=status.HTTP_204_NO_CONTENT)
@limiter.limit(RateLimits.FILE_DELETE)
async def delete_file(
    request: Request,
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete file

    Removes file from storage and updates user storage quota
    """
    file = db.query(File).join(Sample).join(Project).filter(
        File.id == file_id,
        Project.user_id == current_user.id
    ).first()

    if not file:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="File not found"
        )

    try:
        # Delete from storage
        await storage_service.delete_file(file.file_path)

        # Update user storage
        if file.file_size:
            current_user.storage_used = max(0, current_user.storage_used - file.file_size)

        # Delete record
        db.delete(file)
        db.commit()

        return None

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error deleting file: {str(e)}"
        )
