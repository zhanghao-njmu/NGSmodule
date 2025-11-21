"""
Pipeline Template management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import Optional, List
from uuid import UUID
from collections import defaultdict

from app.core.database import get_db
from app.core.deps import get_current_user, get_current_admin
from app.models.user import User
from app.models.pipeline_template import PipelineTemplate
from app.models.project import Project
from app.models.task import PipelineTask
from app.schemas.pipeline import (
    PipelineTemplateCreate,
    PipelineTemplateUpdate,
    PipelineTemplateResponse,
    PipelineTemplateListResponse,
    PipelineExecuteRequest,
    PipelineTemplateCategory,
)
from app.schemas.task import TaskResponse
from app.schemas.common import MessageResponse
from app.workers.pipeline_tasks import run_ngs_pipeline
from datetime import datetime

router = APIRouter()


@router.get("", response_model=PipelineTemplateListResponse)
async def list_pipeline_templates(
    category: Optional[str] = Query(None, description="Filter by category"),
    is_active: Optional[bool] = Query(None, description="Filter by active status"),
    search: Optional[str] = Query(None, description="Search in name and description"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    List all available pipeline templates
    """
    query = db.query(PipelineTemplate)

    # Apply filters
    if category:
        query = query.filter(PipelineTemplate.category == category)

    if is_active is not None:
        query = query.filter(PipelineTemplate.is_active == is_active)

    if search:
        search_pattern = f"%{search}%"
        query = query.filter(
            (PipelineTemplate.name.ilike(search_pattern)) |
            (PipelineTemplate.display_name.ilike(search_pattern)) |
            (PipelineTemplate.description.ilike(search_pattern))
        )

    # Order by sort_order and name
    templates = query.order_by(
        PipelineTemplate.sort_order,
        PipelineTemplate.display_name
    ).all()

    return PipelineTemplateListResponse(
        total=len(templates),
        items=templates
    )


@router.get("/categories", response_model=List[PipelineTemplateCategory])
async def get_pipeline_categories(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get all pipeline categories with template counts
    """
    templates = db.query(PipelineTemplate).filter(
        PipelineTemplate.is_active == True
    ).all()

    # Group by category
    categories = defaultdict(lambda: {"count": 0, "templates": []})
    for template in templates:
        categories[template.category]["count"] += 1
        categories[template.category]["templates"].append(template.display_name)

    return [
        PipelineTemplateCategory(
            category=cat,
            count=data["count"],
            templates=data["templates"]
        )
        for cat, data in sorted(categories.items())
    ]


@router.get("/{template_id}", response_model=PipelineTemplateResponse)
async def get_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get pipeline template details by ID
    """
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == template_id
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found"
        )

    return template


@router.post("", response_model=PipelineTemplateResponse, status_code=status.HTTP_201_CREATED)
async def create_pipeline_template(
    template_data: PipelineTemplateCreate,
    current_user: User = Depends(get_current_admin),  # Only admin can create templates
    db: Session = Depends(get_db)
):
    """
    Create a new custom pipeline template (Admin only)
    """
    # Check if template name already exists
    existing = db.query(PipelineTemplate).filter(
        PipelineTemplate.name == template_data.name
    ).first()

    if existing:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Pipeline template with name '{template_data.name}' already exists"
        )

    # Create new template
    template = PipelineTemplate(
        **template_data.model_dump(),
        is_active=True,
        is_builtin=False  # Custom templates
    )

    db.add(template)
    db.commit()
    db.refresh(template)

    return template


@router.put("/{template_id}", response_model=PipelineTemplateResponse)
async def update_pipeline_template(
    template_id: UUID,
    template_data: PipelineTemplateUpdate,
    current_user: User = Depends(get_current_admin),  # Only admin can update
    db: Session = Depends(get_db)
):
    """
    Update pipeline template (Admin only)

    Note: Built-in templates have restrictions on what can be updated
    """
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == template_id
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found"
        )

    # Update fields
    update_data = template_data.model_dump(exclude_unset=True)

    # Prevent modifying critical fields for built-in templates
    if template.is_builtin:
        restricted_fields = ['name', 'script_name']
        for field in restricted_fields:
            if field in update_data:
                del update_data[field]

    for field, value in update_data.items():
        setattr(template, field, value)

    db.commit()
    db.refresh(template)

    return template


@router.delete("/{template_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_admin),  # Only admin can delete
    db: Session = Depends(get_db)
):
    """
    Delete pipeline template (Admin only)

    Note: Built-in templates cannot be deleted
    """
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == template_id
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found"
        )

    if template.is_builtin:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot delete built-in pipeline templates"
        )

    db.delete(template)
    db.commit()

    return None


@router.post("/execute", response_model=TaskResponse, status_code=status.HTTP_201_CREATED)
async def execute_pipeline(
    execute_data: PipelineExecuteRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Execute a pipeline with specified parameters

    Creates a new task and starts pipeline execution using Celery
    """
    # Verify template exists
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == execute_data.template_id,
        PipelineTemplate.is_active == True
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found or inactive"
        )

    # Verify project belongs to user
    project = db.query(Project).filter(
        Project.id == execute_data.project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Merge default parameters with user-provided parameters
    merged_params = {**template.default_params, **execute_data.parameters}

    # Add template information to config
    config = {
        "template_id": str(execute_data.template_id),
        "template_name": template.name,
        "script_name": template.script_name,
        "script_path": template.script_path,
        "sample_ids": [str(sid) for sid in execute_data.sample_ids],
        "parameters": merged_params
    }

    # Create task
    task = PipelineTask(
        project_id=execute_data.project_id,
        task_name=execute_data.task_name,
        task_type=template.category,
        config=config,
        status="pending"
    )

    db.add(task)
    db.commit()
    db.refresh(task)

    # Submit to Celery for execution
    try:
        celery_task = run_ngs_pipeline.delay(
            task_id=str(task.id),
            pipeline_script=template.script_name,
            config=config
        )

        task.celery_task_id = celery_task.id
        task.status = "pending"
        task.started_at = datetime.utcnow()
        db.commit()
        db.refresh(task)

        return task

    except Exception as e:
        task.status = "failed"
        task.error_message = f"Failed to start pipeline: {str(e)}"
        db.commit()

        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error starting pipeline: {str(e)}"
        )


@router.post("/{template_id}/toggle", response_model=MessageResponse)
async def toggle_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Toggle pipeline template active status (Admin only)
    """
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == template_id
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found"
        )

    template.is_active = not template.is_active
    db.commit()

    status_text = "activated" if template.is_active else "deactivated"
    return MessageResponse(message=f"Pipeline template {status_text} successfully")
