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
from app.models.sample import Sample
from app.schemas.pipeline import (
    PipelineTemplateCreate,
    PipelineTemplateUpdate,
    PipelineTemplateResponse,
    PipelineTemplateListResponse,
    PipelineExecuteRequest,
    PipelineTemplateCategory,
    PipelineBatchExecuteRequest,
    PipelineBatchExecuteResponse,
    ParameterRecommendationResponse,
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


@router.post("/batch-execute", response_model=PipelineBatchExecuteResponse, status_code=status.HTTP_201_CREATED)
async def batch_execute_pipeline(
    execute_data: PipelineBatchExecuteRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Execute a pipeline on multiple samples in batch

    Creates one task per sample for parallel processing
    """
    # Verify template exists and is active
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

    # Verify all samples belong to the project
    samples = db.query(Sample).filter(
        Sample.id.in_(execute_data.sample_ids),
        Sample.project_id == execute_data.project_id
    ).all()

    if len(samples) != len(execute_data.sample_ids):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Some samples not found or don't belong to the project"
        )

    # Merge default parameters with user-provided parameters
    merged_params = {**template.default_params, **execute_data.parameters}

    created_tasks = []
    failed_samples = []

    # Create and execute task for each sample
    for sample in samples:
        try:
            # Create task name with sample identifier
            task_name = f"{execute_data.task_name_prefix} - {sample.name}"

            # Add template and sample information to config
            config = {
                "template_id": str(execute_data.template_id),
                "template_name": template.name,
                "script_name": template.script_name,
                "script_path": template.script_path,
                "sample_ids": [str(sample.id)],
                "sample_name": sample.name,
                "parameters": merged_params
            }

            # Create task
            task = PipelineTask(
                project_id=execute_data.project_id,
                task_name=task_name,
                task_type=template.category,
                config=config,
                status="pending"
            )

            db.add(task)
            db.flush()  # Get task ID without committing

            # Submit to Celery for execution
            celery_task = run_ngs_pipeline.delay(
                task_id=str(task.id),
                pipeline_script=template.script_name,
                config=config
            )

            task.celery_task_id = celery_task.id
            task.status = "pending"
            task.started_at = datetime.utcnow()

            created_tasks.append(task.id)

        except Exception as e:
            failed_samples.append({
                "sample_id": str(sample.id),
                "sample_name": sample.name,
                "error": str(e)
            })

    # Commit all successful tasks
    db.commit()

    return PipelineBatchExecuteResponse(
        total_tasks=len(created_tasks),
        created_tasks=created_tasks,
        failed_samples=failed_samples
    )


@router.get("/{template_id}/recommend-parameters", response_model=ParameterRecommendationResponse)
async def recommend_parameters(
    template_id: UUID,
    project_id: Optional[UUID] = Query(None, description="Filter by project for personalized recommendations"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get AI-powered parameter recommendations based on historical successful tasks

    Analyzes completed tasks using the same template to recommend optimal parameters
    """
    # Verify template exists
    template = db.query(PipelineTemplate).filter(
        PipelineTemplate.id == template_id,
        PipelineTemplate.is_active == True
    ).first()

    if not template:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Pipeline template not found or inactive"
        )

    # Query successful tasks for this template
    query = db.query(PipelineTask).filter(
        PipelineTask.config['template_id'].astext == str(template_id),
        PipelineTask.status == 'completed'
    )

    # If project_id provided, prioritize tasks from the same project
    if project_id:
        # Verify user has access to the project
        project = db.query(Project).filter(
            Project.id == project_id,
            Project.user_id == current_user.id
        ).first()

        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Project not found"
            )

        # Get project-specific tasks first
        project_tasks = query.filter(PipelineTask.project_id == project_id).limit(20).all()
        if len(project_tasks) < 5:
            # If not enough project-specific tasks, include all user's tasks
            query = query.filter(PipelineTask.project_id.in_(
                db.query(Project.id).filter(Project.user_id == current_user.id)
            ))
            all_tasks = query.limit(50).all()
        else:
            all_tasks = project_tasks
    else:
        # Get tasks from all user's projects
        query = query.filter(PipelineTask.project_id.in_(
            db.query(Project.id).filter(Project.user_id == current_user.id)
        ))
        all_tasks = query.limit(50).all()

    # If no historical tasks found, return template defaults
    if not all_tasks:
        return ParameterRecommendationResponse(
            recommended_params=template.default_params,
            confidence_score=0.5,
            based_on_tasks=0,
            explanation="No historical tasks found. Using template default parameters."
        )

    # Analyze parameters from successful tasks
    param_stats = defaultdict(lambda: defaultdict(int))
    total_tasks = len(all_tasks)

    for task in all_tasks:
        if 'parameters' in task.config:
            params = task.config['parameters']
            for key, value in params.items():
                # Count frequency of each parameter value
                param_stats[key][str(value)] += 1

    # Build recommended parameters
    recommended_params = {}
    for key, value_counts in param_stats.items():
        # Get most common value
        most_common_value = max(value_counts.items(), key=lambda x: x[1])
        value_str, count = most_common_value

        # Try to convert back to appropriate type
        try:
            if value_str.lower() in ['true', 'false']:
                recommended_params[key] = value_str.lower() == 'true'
            elif '.' in value_str:
                recommended_params[key] = float(value_str)
            elif value_str.isdigit():
                recommended_params[key] = int(value_str)
            else:
                recommended_params[key] = value_str
        except:
            recommended_params[key] = value_str

    # Merge with template defaults for missing parameters
    final_params = {**template.default_params, **recommended_params}

    # Calculate confidence score based on task count and parameter consistency
    if total_tasks >= 20:
        base_confidence = 0.9
    elif total_tasks >= 10:
        base_confidence = 0.8
    elif total_tasks >= 5:
        base_confidence = 0.7
    else:
        base_confidence = 0.6

    # Adjust confidence based on parameter consistency
    avg_consistency = sum(
        max(counts.values()) / sum(counts.values())
        for counts in param_stats.values()
    ) / len(param_stats) if param_stats else 0.5

    confidence_score = min(base_confidence * avg_consistency, 1.0)

    # Generate explanation
    if project_id and all_tasks[0].project_id == project_id:
        explanation = f"Recommendations based on {total_tasks} successful tasks from this project. "
    else:
        explanation = f"Recommendations based on {total_tasks} successful tasks from your projects. "

    explanation += f"These parameters were used in {int(avg_consistency * 100)}% of successful runs."

    return ParameterRecommendationResponse(
        recommended_params=final_params,
        confidence_score=round(confidence_score, 2),
        based_on_tasks=total_tasks,
        explanation=explanation
    )
