"""
Pipeline Task management API endpoints (Refactored with Service Layer)

This is the refactored version using TaskService for business logic.
Routers are now thin HTTP adapters that delegate to the service layer.
"""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Query, status
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.schemas.common import MessageResponse
from app.schemas.task import (
    TaskCreate,
    TaskExecuteRequest,
    TaskListResponse,
    TaskLogResponse,
    TaskResponse,
    TaskStats,
    TaskUpdate,
)
from app.services.task_service import TaskService

router = APIRouter()


# ============= DEPENDENCY INJECTION =============


def get_task_service(db: Session = Depends(get_db)) -> TaskService:
    """Dependency to get TaskService instance"""
    return TaskService(db)


# ============= ENDPOINTS =============


@router.get("", response_model=TaskListResponse)
async def list_tasks(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    status: Optional[str] = Query(None, description="Filter by status"),
    task_type: Optional[str] = Query(None, description="Filter by task type"),
    current_user: User = Depends(get_current_user),
    service: TaskService = Depends(get_task_service),
):
    """
    List tasks with optional filters
    """
    tasks = service.list_tasks(
        user_id=current_user.id, project_id=project_id, status_filter=status, task_type=task_type
    )

    return TaskListResponse(total=len(tasks), items=tasks)


@router.get("/stats", response_model=TaskStats)
async def get_task_stats(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    current_user: User = Depends(get_current_user),
    service: TaskService = Depends(get_task_service),
):
    """
    Get task statistics
    """
    return service.get_stats(user_id=current_user.id, project_id=project_id)


@router.post("", response_model=TaskResponse, status_code=status.HTTP_201_CREATED)
async def create_task(
    task_data: TaskCreate,
    current_user: User = Depends(get_current_user),
    service: TaskService = Depends(get_task_service),
):
    """
    Create a new pipeline task

    - **task_name**: Name of the task
    - **task_type**: Type of task (e.g., 'RNA-seq', 'DNA-seq', 'scRNA-seq')
    - **project_id**: Project ID to associate task with
    - **config**: Task configuration parameters
    """
    return service.create(current_user.id, task_data)


@router.get("/{task_id}", response_model=TaskResponse)
async def get_task(
    task_id: UUID, current_user: User = Depends(get_current_user), service: TaskService = Depends(get_task_service)
):
    """
    Get task details by ID
    """
    return service.get_by_id_or_raise(task_id, current_user.id)


@router.put("/{task_id}", response_model=TaskResponse)
async def update_task(
    task_id: UUID,
    task_data: TaskUpdate,
    current_user: User = Depends(get_current_user),
    service: TaskService = Depends(get_task_service),
):
    """
    Update task details

    Can update task name, type, status, progress, error message, or configuration
    """
    return service.update(task_id, current_user.id, task_data)


@router.post("/{task_id}/execute", response_model=MessageResponse)
async def execute_task(
    task_id: UUID,
    execute_data: TaskExecuteRequest,
    current_user: User = Depends(get_current_user),
    service: TaskService = Depends(get_task_service),
):
    """
    Execute a pipeline task

    Triggers Celery worker to run the NGS pipeline script with the provided configuration

    - **pipeline_script**: Name of the pipeline script to execute (e.g., 'RNAseq_pipeline.sh')
    - **config**: Pipeline configuration parameters (e.g., thread count, quality thresholds)
    """
    message, celery_task_id = service.execute(task_id, current_user.id, execute_data)
    return MessageResponse(message=message)


@router.post("/{task_id}/cancel", response_model=MessageResponse)
async def cancel_task(
    task_id: UUID, current_user: User = Depends(get_current_user), service: TaskService = Depends(get_task_service)
):
    """
    Cancel a running task

    Attempts to revoke the Celery task and marks it as cancelled
    """
    message = service.cancel(task_id, current_user.id)
    return MessageResponse(message=message)


@router.get("/{task_id}/logs", response_model=TaskLogResponse)
async def get_task_logs(
    task_id: UUID, current_user: User = Depends(get_current_user), service: TaskService = Depends(get_task_service)
):
    """
    Get task execution logs

    Returns the content of the task's log file if available
    """
    log_content, log_file_path = service.get_logs(task_id, current_user.id)

    return TaskLogResponse(task_id=task_id, log_content=log_content, log_file_path=log_file_path)


@router.delete("/{task_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_task(
    task_id: UUID, current_user: User = Depends(get_current_user), service: TaskService = Depends(get_task_service)
):
    """
    Delete a task

    Only tasks with status 'pending', 'completed', 'failed', or 'cancelled' can be deleted.
    Running tasks must be cancelled first.
    """
    service.delete(task_id, current_user.id)
    return None
