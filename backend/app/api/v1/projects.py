"""
Project management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from sqlalchemy import func
from typing import List, Optional
from uuid import UUID

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.schemas.project import (
    ProjectCreate,
    ProjectUpdate,
    ProjectResponse,
    ProjectListResponse,
    ProjectStats,
)
from app.schemas.common import MessageResponse

router = APIRouter()


@router.get("/stats", response_model=ProjectStats)
async def get_project_stats(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get project statistics for current user
    """
    # Count projects
    total_projects = db.query(Project).filter(Project.user_id == current_user.id).count()
    active_projects = db.query(Project).filter(
        Project.user_id == current_user.id,
        Project.status == "active"
    ).count()
    archived_projects = db.query(Project).filter(
        Project.user_id == current_user.id,
        Project.status == "archived"
    ).count()

    # Count samples
    total_samples = db.query(Sample).join(Project).filter(
        Project.user_id == current_user.id
    ).count()

    # Count tasks
    total_tasks = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id
    ).count()
    running_tasks = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id,
        PipelineTask.status == "running"
    ).count()
    completed_tasks = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id,
        PipelineTask.status == "completed"
    ).count()
    failed_tasks = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id,
        PipelineTask.status == "failed"
    ).count()

    return ProjectStats(
        total_projects=total_projects,
        active_projects=active_projects,
        archived_projects=archived_projects,
        total_samples=total_samples,
        total_tasks=total_tasks,
        running_tasks=running_tasks,
        completed_tasks=completed_tasks,
        failed_tasks=failed_tasks,
    )


@router.get("", response_model=ProjectListResponse)
async def list_projects(
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    status: Optional[str] = Query(None, description="Filter by status"),
    project_type: Optional[str] = Query(None, description="Filter by type"),
    search: Optional[str] = Query(None, description="Search by name"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    List projects for current user with pagination and filters
    """
    query = db.query(Project).filter(Project.user_id == current_user.id)

    # Apply filters
    if status:
        query = query.filter(Project.status == status)
    if project_type:
        query = query.filter(Project.project_type == project_type)
    if search:
        query = query.filter(Project.name.ilike(f"%{search}%"))

    # Get total count
    total = query.count()

    # Get paginated results
    projects = query.order_by(Project.created_at.desc()).offset(skip).limit(limit).all()

    # Add computed fields
    for project in projects:
        project.sample_count = db.query(Sample).filter(Sample.project_id == project.id).count()
        project.task_count = db.query(PipelineTask).filter(PipelineTask.project_id == project.id).count()

    return ProjectListResponse(
        total=total,
        items=projects,
        page=skip // limit + 1,
        page_size=limit
    )


@router.post("", response_model=ProjectResponse, status_code=status.HTTP_201_CREATED)
async def create_project(
    project_data: ProjectCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Create a new project

    - **name**: Project name (required, unique per user)
    - **description**: Optional project description
    - **project_type**: Type of analysis (rna-seq, dna-seq, etc.)
    - **config**: Optional project configuration
    """
    # Check if project name already exists for this user
    existing = db.query(Project).filter(
        Project.user_id == current_user.id,
        Project.name == project_data.name
    ).first()

    if existing:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Project '{project_data.name}' already exists"
        )

    # Create project
    project = Project(
        user_id=current_user.id,
        name=project_data.name,
        description=project_data.description,
        project_type=project_data.project_type,
        config=project_data.config,
        status="active"
    )

    db.add(project)
    db.commit()
    db.refresh(project)

    # Add computed fields
    project.sample_count = 0
    project.task_count = 0

    return project


@router.get("/{project_id}", response_model=ProjectResponse)
async def get_project(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get project by ID
    """
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Add computed fields
    project.sample_count = db.query(Sample).filter(Sample.project_id == project.id).count()
    project.task_count = db.query(PipelineTask).filter(PipelineTask.project_id == project.id).count()

    return project


@router.put("/{project_id}", response_model=ProjectResponse)
async def update_project(
    project_id: UUID,
    project_update: ProjectUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Update project

    - **name**: New project name
    - **description**: New description
    - **project_type**: New project type
    - **status**: New status (active, archived)
    - **config**: Updated configuration
    """
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Check if new name conflicts
    if project_update.name and project_update.name != project.name:
        existing = db.query(Project).filter(
            Project.user_id == current_user.id,
            Project.name == project_update.name,
            Project.id != project_id
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Project '{project_update.name}' already exists"
            )

    # Update fields
    update_data = project_update.dict(exclude_unset=True)
    for field, value in update_data.items():
        setattr(project, field, value)

    db.commit()
    db.refresh(project)

    # Add computed fields
    project.sample_count = db.query(Sample).filter(Sample.project_id == project.id).count()
    project.task_count = db.query(PipelineTask).filter(PipelineTask.project_id == project.id).count()

    return project


@router.delete("/{project_id}", response_model=MessageResponse)
async def delete_project(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete project (soft delete by setting status to 'deleted')
    """
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Soft delete
    project.status = "deleted"
    db.commit()

    return MessageResponse(
        message="Project deleted successfully",
        detail=f"Project '{project.name}' has been marked as deleted"
    )


@router.post("/{project_id}/archive", response_model=ProjectResponse)
async def archive_project(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Archive project
    """
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    project.status = "archived"
    db.commit()
    db.refresh(project)

    return project


@router.post("/{project_id}/restore", response_model=ProjectResponse)
async def restore_project(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Restore archived or deleted project
    """
    project = db.query(Project).filter(
        Project.id == project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    project.status = "active"
    db.commit()
    db.refresh(project)

    return project
