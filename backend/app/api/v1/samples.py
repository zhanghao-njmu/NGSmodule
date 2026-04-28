"""
Sample management API endpoints (Refactored with Service Layer)

This is the refactored version using SampleService for business logic.
Routers are now thin HTTP adapters that delegate to the service layer.
"""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, File, Query, Request, UploadFile, status
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_user
from app.core.rate_limit import RateLimits, limiter
from app.models.user import User
from app.schemas.common import MessageResponse
from app.schemas.sample import (
    SampleBatchCreate,
    SampleCreate,
    SampleListResponse,
    SampleResponse,
    SampleUpdate,
)
from app.services.sample_service import SampleService

router = APIRouter()


# ============= DEPENDENCY INJECTION =============


def get_sample_service(db: Session = Depends(get_db)) -> SampleService:
    """Dependency to get SampleService instance"""
    return SampleService(db)


# ============= ENDPOINTS =============


@router.get("", response_model=SampleListResponse)
async def list_samples(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    group_name: Optional[str] = Query(None, description="Filter by group"),
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=500),
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    List samples with optional filters
    """
    samples, total = service.list_samples(
        user_id=current_user.id, skip=skip, limit=limit, project_id=project_id, group_name=group_name
    )

    return SampleListResponse(total=total, items=samples, page=skip // limit + 1, page_size=limit)


@router.post("", response_model=SampleResponse, status_code=status.HTTP_201_CREATED)
async def create_sample(
    sample_data: SampleCreate,
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Create a new sample

    - **project_id**: Project ID (required)
    - **sample_id**: Sample identifier (required, unique per project)
    - **run_id**: Sequencing run ID
    - **group_name**: Sample group (e.g., control, treatment)
    - **layout**: Sequencing layout (PE or SE)
    - **batch_id**: Batch identifier
    - **metadata**: Additional sample metadata
    """
    return service.create(current_user.id, sample_data)


@router.post("/batch", response_model=MessageResponse, status_code=status.HTTP_201_CREATED)
async def create_samples_batch(
    batch_data: SampleBatchCreate,
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Create multiple samples at once
    """
    count = service.create_batch(
        user_id=current_user.id, project_id=batch_data.project_id, samples_data=batch_data.samples
    )

    return MessageResponse(message=f"Successfully created {count} samples", detail="Added to project")


@router.post("/import-csv", response_model=MessageResponse, status_code=status.HTTP_201_CREATED)
@limiter.limit(RateLimits.BATCH_IMPORT)
async def import_samples_from_csv(
    request: Request,
    project_id: UUID,
    file: UploadFile = File(..., description="CSV file with sample information"),
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Import samples from CSV file

    CSV format (with headers):
    ```
    sample_id,run_id,group_name,layout,batch_id
    Sample1,Run001,Control,PE,Batch1
    Sample2,Run002,Treatment,PE,Batch1
    ```
    """
    # Read CSV content
    content = await file.read()

    count, project_name = service.import_from_csv(user_id=current_user.id, project_id=project_id, csv_content=content)

    return MessageResponse(
        message=f"Successfully imported {count} samples", detail=f"Added to project '{project_name}'"
    )


@router.get("/{sample_id}", response_model=SampleResponse)
async def get_sample(
    sample_id: UUID,
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Get sample by ID
    """
    return service.get_by_id_or_raise(sample_id, current_user.id)


@router.put("/{sample_id}", response_model=SampleResponse)
async def update_sample(
    sample_id: UUID,
    sample_update: SampleUpdate,
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Update sample information
    """
    return service.update(sample_id, current_user.id, sample_update)


@router.delete("/{sample_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_sample(
    sample_id: UUID,
    current_user: User = Depends(get_current_user),
    service: SampleService = Depends(get_sample_service),
):
    """
    Delete sample and all associated files
    """
    service.delete(sample_id, current_user.id)
    return None
