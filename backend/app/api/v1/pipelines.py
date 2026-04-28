"""
Pipeline Template management API endpoints (Refactored with Service Layer)

This is the refactored version using PipelineService for business logic.
Routers are now thin HTTP adapters that delegate to the service layer.
"""

from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Query, status
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_admin, get_current_user
from app.models.user import User
from app.schemas.common import MessageResponse
from app.schemas.pipeline import (
    ParameterRecommendationResponse,
    PipelineBatchExecuteRequest,
    PipelineBatchExecuteResponse,
    PipelineExecuteRequest,
    PipelineTemplateCategory,
    PipelineTemplateCreate,
    PipelineTemplateListResponse,
    PipelineTemplateResponse,
    PipelineTemplateUpdate,
)
from app.schemas.task import TaskResponse
from app.services.pipeline_service import PipelineService

router = APIRouter()


# ============= DEPENDENCY INJECTION =============


def get_pipeline_service(db: Session = Depends(get_db)) -> PipelineService:
    """Dependency to get PipelineService instance"""
    return PipelineService(db)


# ============= ENDPOINTS =============


@router.get("", response_model=PipelineTemplateListResponse)
async def list_pipeline_templates(
    category: Optional[str] = Query(None, description="Filter by category"),
    is_active: Optional[bool] = Query(None, description="Filter by active status"),
    search: Optional[str] = Query(None, description="Search in name and description"),
    current_user: User = Depends(get_current_user),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    List all available pipeline templates
    """
    templates = service.list_templates(category=category, is_active=is_active, search=search)

    return PipelineTemplateListResponse(total=len(templates), items=templates)


@router.get("/categories", response_model=List[PipelineTemplateCategory])
async def get_pipeline_categories(
    current_user: User = Depends(get_current_user), service: PipelineService = Depends(get_pipeline_service)
):
    """
    Get all pipeline categories with template counts
    """
    return service.get_categories()


@router.get("/{template_id}", response_model=PipelineTemplateResponse)
async def get_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_user),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Get pipeline template details by ID
    """
    return service.get_by_id_or_raise(template_id)


@router.post("", response_model=PipelineTemplateResponse, status_code=status.HTTP_201_CREATED)
async def create_pipeline_template(
    template_data: PipelineTemplateCreate,
    current_user: User = Depends(get_current_admin),  # Only admin can create templates
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Create a new custom pipeline template (Admin only)
    """
    return service.create(template_data)


@router.put("/{template_id}", response_model=PipelineTemplateResponse)
async def update_pipeline_template(
    template_id: UUID,
    template_data: PipelineTemplateUpdate,
    current_user: User = Depends(get_current_admin),  # Only admin can update
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Update pipeline template (Admin only)

    Note: Built-in templates have restrictions on what can be updated
    """
    return service.update(template_id, template_data)


@router.delete("/{template_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_admin),  # Only admin can delete
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Delete pipeline template (Admin only)

    Note: Built-in templates cannot be deleted
    """
    service.delete(template_id)
    return None


@router.post("/{template_id}/toggle", response_model=MessageResponse)
async def toggle_pipeline_template(
    template_id: UUID,
    current_user: User = Depends(get_current_admin),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Toggle pipeline template active status (Admin only)
    """
    template, message = service.toggle_active(template_id)
    return MessageResponse(message=message)


@router.post("/execute", response_model=TaskResponse, status_code=status.HTTP_201_CREATED)
async def execute_pipeline(
    execute_data: PipelineExecuteRequest,
    current_user: User = Depends(get_current_user),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Execute a pipeline with specified parameters

    Creates a new task and starts pipeline execution using Celery
    """
    return service.execute(current_user.id, execute_data)


@router.post("/batch-execute", response_model=PipelineBatchExecuteResponse, status_code=status.HTTP_201_CREATED)
async def batch_execute_pipeline(
    execute_data: PipelineBatchExecuteRequest,
    current_user: User = Depends(get_current_user),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Execute a pipeline on multiple samples in batch

    Creates one task per sample for parallel processing
    """
    return service.batch_execute(current_user.id, execute_data)


@router.get("/{template_id}/recommend-parameters", response_model=ParameterRecommendationResponse)
async def recommend_parameters(
    template_id: UUID,
    project_id: Optional[UUID] = Query(None, description="Filter by project for personalized recommendations"),
    current_user: User = Depends(get_current_user),
    service: PipelineService = Depends(get_pipeline_service),
):
    """
    Get AI-powered parameter recommendations based on historical successful tasks

    Analyzes completed tasks using the same template to recommend optimal parameters
    """
    return service.recommend_parameters(template_id=template_id, user_id=current_user.id, project_id=project_id)
