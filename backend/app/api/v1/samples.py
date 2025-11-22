"""
Sample management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query, UploadFile, File
from sqlalchemy.orm import Session
from sqlalchemy import func
from typing import List, Optional
from uuid import UUID
import csv
import io

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.file import File as FileModel
from app.schemas.sample import (
    SampleCreate,
    SampleBatchCreate,
    SampleUpdate,
    SampleResponse,
    SampleListResponse,
)
from app.schemas.common import MessageResponse

router = APIRouter()


@router.get("", response_model=SampleListResponse)
async def list_samples(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    group_name: Optional[str] = Query(None, description="Filter by group"),
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=500),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    List samples with optional filters
    """
    query = db.query(Sample).join(Project).filter(Project.user_id == current_user.id)

    # Apply filters
    if project_id:
        # Verify project belongs to user
        project = db.query(Project).filter(
            Project.id == project_id,
            Project.user_id == current_user.id
        ).first()
        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Project not found"
            )
        query = query.filter(Sample.project_id == project_id)

    if group_name:
        query = query.filter(Sample.group_name == group_name)

    # Get total count
    total = query.count()

    # Get paginated results
    samples = query.order_by(Sample.created_at.desc()).offset(skip).limit(limit).all()

    # Add computed fields - Optimized to avoid N+1 queries
    if samples:
        sample_ids = [s.id for s in samples]

        # Get file counts in one query
        file_counts = dict(
            db.query(FileModel.sample_id, func.count(FileModel.id))
            .filter(FileModel.sample_id.in_(sample_ids))
            .group_by(FileModel.sample_id)
            .all()
        )

        # Assign counts to samples
        for sample in samples:
            sample.file_count = file_counts.get(sample.id, 0)

    return SampleListResponse(
        total=total,
        items=samples,
        page=skip // limit + 1,
        page_size=limit
    )


@router.post("", response_model=SampleResponse, status_code=status.HTTP_201_CREATED)
async def create_sample(
    sample_data: SampleCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
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
    # Verify project belongs to user
    project = db.query(Project).filter(
        Project.id == sample_data.project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Check if sample_id already exists in this project
    existing = db.query(Sample).filter(
        Sample.project_id == sample_data.project_id,
        Sample.sample_id == sample_data.sample_id
    ).first()

    if existing:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Sample '{sample_data.sample_id}' already exists in this project"
        )

    # Create sample
    sample = Sample(
        project_id=sample_data.project_id,
        sample_id=sample_data.sample_id,
        run_id=sample_data.run_id,
        group_name=sample_data.group_name,
        layout=sample_data.layout,
        batch_id=sample_data.batch_id,
        metadata=sample_data.metadata or {}
    )

    db.add(sample)
    db.commit()
    db.refresh(sample)

    # Add computed fields
    sample.file_count = 0

    return sample


@router.post("/batch", response_model=MessageResponse, status_code=status.HTTP_201_CREATED)
async def create_samples_batch(
    batch_data: SampleBatchCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Create multiple samples at once
    """
    # Verify project belongs to user
    project = db.query(Project).filter(
        Project.id == batch_data.project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Get existing sample IDs
    existing_samples = db.query(Sample.sample_id).filter(
        Sample.project_id == batch_data.project_id
    ).all()
    existing_ids = {s.sample_id for s in existing_samples}

    # Check for duplicates
    new_ids = [s.sample_id for s in batch_data.samples]
    duplicates = [sid for sid in new_ids if sid in existing_ids]

    if duplicates:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Duplicate sample IDs found: {', '.join(duplicates)}"
        )

    # Create samples
    samples_to_add = []
    for sample_data in batch_data.samples:
        sample = Sample(
            project_id=batch_data.project_id,
            sample_id=sample_data.sample_id,
            run_id=sample_data.run_id,
            group_name=sample_data.group_name,
            layout=sample_data.layout,
            batch_id=sample_data.batch_id,
            metadata={}
        )
        samples_to_add.append(sample)

    db.add_all(samples_to_add)
    db.commit()

    return MessageResponse(
        message=f"Successfully created {len(samples_to_add)} samples",
        detail=f"Added to project '{project.name}'"
    )


@router.post("/import-csv", response_model=MessageResponse, status_code=status.HTTP_201_CREATED)
async def import_samples_from_csv(
    project_id: UUID,
    file: UploadFile = File(..., description="CSV file with sample information"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
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
    # Verify project belongs to user
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Read CSV content
    try:
        content = await file.read()
        decoded = content.decode('utf-8')
        csv_reader = csv.DictReader(io.StringIO(decoded))

        # Validate headers
        required_headers = {'sample_id'}
        if not required_headers.issubset(csv_reader.fieldnames or []):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"CSV must contain at least: {', '.join(required_headers)}"
            )

        # Get existing sample IDs
        existing_samples = db.query(Sample.sample_id).filter(
            Sample.project_id == project_id
        ).all()
        existing_ids = {s.sample_id for s in existing_samples}

        # Parse and create samples
        samples_to_add = []
        for row in csv_reader:
            sample_id = row.get('sample_id', '').strip()
            if not sample_id:
                continue

            if sample_id in existing_ids:
                continue  # Skip duplicates

            sample = Sample(
                project_id=project_id,
                sample_id=sample_id,
                run_id=row.get('run_id', '').strip() or None,
                group_name=row.get('group_name', '').strip() or None,
                layout=row.get('layout', '').strip() or None,
                batch_id=row.get('batch_id', '').strip() or None,
                metadata={}
            )
            samples_to_add.append(sample)
            existing_ids.add(sample_id)

        if not samples_to_add:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="No valid samples found in CSV or all samples already exist"
            )

        db.add_all(samples_to_add)
        db.commit()

        return MessageResponse(
            message=f"Successfully imported {len(samples_to_add)} samples",
            detail=f"Added to project '{project.name}'"
        )

    except UnicodeDecodeError:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid CSV file encoding. Please use UTF-8 encoding."
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Error parsing CSV: {str(e)}"
        )


@router.get("/{sample_id}", response_model=SampleResponse)
async def get_sample(
    sample_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get sample by ID
    """
    sample = db.query(Sample).join(Project).filter(
        Sample.id == sample_id,
        Project.user_id == current_user.id
    ).first()

    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Sample not found"
        )

    # Add computed fields
    sample.file_count = db.query(FileModel).filter(FileModel.sample_id == sample.id).count()

    return sample


@router.put("/{sample_id}", response_model=SampleResponse)
async def update_sample(
    sample_id: UUID,
    sample_update: SampleUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Update sample information
    """
    sample = db.query(Sample).join(Project).filter(
        Sample.id == sample_id,
        Project.user_id == current_user.id
    ).first()

    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Sample not found"
        )

    # Check if new sample_id conflicts
    if sample_update.sample_id and sample_update.sample_id != sample.sample_id:
        existing = db.query(Sample).filter(
            Sample.project_id == sample.project_id,
            Sample.sample_id == sample_update.sample_id,
            Sample.id != sample_id
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Sample '{sample_update.sample_id}' already exists in this project"
            )

    # Update fields
    update_data = sample_update.dict(exclude_unset=True)
    for field, value in update_data.items():
        setattr(sample, field, value)

    db.commit()
    db.refresh(sample)

    # Add computed fields
    sample.file_count = db.query(FileModel).filter(FileModel.sample_id == sample.id).count()

    return sample


@router.delete("/{sample_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_sample(
    sample_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete sample and all associated files
    """
    sample = db.query(Sample).join(Project).filter(
        Sample.id == sample_id,
        Project.user_id == current_user.id
    ).first()

    if not sample:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Sample not found"
        )

    # Delete sample (cascade will delete associated files)
    db.delete(sample)
    db.commit()

    return None
