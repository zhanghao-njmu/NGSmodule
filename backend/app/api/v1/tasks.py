"""
Pipeline Task management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from typing import Optional
from uuid import UUID
from datetime import datetime
from pathlib import Path

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.models.project import Project
from app.models.task import PipelineTask
from app.schemas.task import (
    TaskCreate,
    TaskUpdate,
    TaskResponse,
    TaskListResponse,
    TaskStats,
    TaskLogResponse,
    TaskExecuteRequest,
)
from app.schemas.common import MessageResponse
from app.workers.pipeline_tasks import run_ngs_pipeline

router = APIRouter()


@router.get("", response_model=TaskListResponse)
async def list_tasks(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    status: Optional[str] = Query(None, description="Filter by status"),
    task_type: Optional[str] = Query(None, description="Filter by task type"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    List tasks with optional filters
    """
    query = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id
    )

    # Apply filters
    if project_id:
        query = query.filter(PipelineTask.project_id == project_id)

    if status:
        query = query.filter(PipelineTask.status == status)

    if task_type:
        query = query.filter(PipelineTask.task_type == task_type)

    tasks = query.order_by(PipelineTask.created_at.desc()).all()

    return TaskListResponse(
        total=len(tasks),
        items=tasks
    )


@router.get("/stats", response_model=TaskStats)
async def get_task_stats(
    project_id: Optional[UUID] = Query(None, description="Filter by project"),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get task statistics
    """
    query = db.query(PipelineTask).join(Project).filter(
        Project.user_id == current_user.id
    )

    if project_id:
        query = query.filter(PipelineTask.project_id == project_id)

    total_tasks = query.count()
    pending_tasks = query.filter(PipelineTask.status == "pending").count()
    running_tasks = query.filter(PipelineTask.status == "running").count()
    completed_tasks = query.filter(PipelineTask.status == "completed").count()
    failed_tasks = query.filter(PipelineTask.status == "failed").count()
    cancelled_tasks = query.filter(PipelineTask.status == "cancelled").count()

    return TaskStats(
        total_tasks=total_tasks,
        pending_tasks=pending_tasks,
        running_tasks=running_tasks,
        completed_tasks=completed_tasks,
        failed_tasks=failed_tasks,
        cancelled_tasks=cancelled_tasks
    )


@router.post("", response_model=TaskResponse, status_code=status.HTTP_201_CREATED)
async def create_task(
    task_data: TaskCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Create a new pipeline task

    - **task_name**: Name of the task
    - **task_type**: Type of task (e.g., 'RNA-seq', 'DNA-seq', 'scRNA-seq')
    - **project_id**: Project ID to associate task with
    - **config**: Task configuration parameters
    """
    # Verify project belongs to user
    project = db.query(Project).filter(
        Project.id == task_data.project_id,
        Project.user_id == current_user.id
    ).first()

    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Project not found"
        )

    # Create task
    task = PipelineTask(
        project_id=task_data.project_id,
        task_name=task_data.task_name,
        task_type=task_data.task_type,
        config=task_data.config,
        status="pending"
    )

    db.add(task)
    db.commit()
    db.refresh(task)

    return task


@router.get("/{task_id}", response_model=TaskResponse)
async def get_task(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get task details by ID
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    return task


@router.put("/{task_id}", response_model=TaskResponse)
async def update_task(
    task_id: UUID,
    task_data: TaskUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Update task details

    Can update task name, type, status, progress, error message, or configuration
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    # Update fields
    update_data = task_data.model_dump(exclude_unset=True)
    for field, value in update_data.items():
        setattr(task, field, value)

    db.commit()
    db.refresh(task)

    return task


@router.post("/{task_id}/execute", response_model=MessageResponse)
async def execute_task(
    task_id: UUID,
    execute_data: TaskExecuteRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Execute a pipeline task

    Triggers Celery worker to run the NGS pipeline script with the provided configuration

    - **pipeline_script**: Name of the pipeline script to execute (e.g., 'RNAseq_pipeline.sh')
    - **config**: Pipeline configuration parameters (e.g., thread count, quality thresholds)
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    # Check if task is already running or completed
    if task.status in ["running", "completed"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Task is already {task.status}"
        )

    # Update task configuration
    task.config = {**task.config, **execute_data.config}
    task.status = "pending"
    task.started_at = datetime.utcnow()
    db.commit()

    # Submit to Celery
    try:
        celery_task = run_ngs_pipeline.delay(
            task_id=str(task.id),
            pipeline_script=execute_data.pipeline_script,
            config=task.config
        )

        task.celery_task_id = celery_task.id
        db.commit()

        return MessageResponse(message=f"Task execution started with Celery ID: {celery_task.id}")

    except Exception as e:
        task.status = "failed"
        task.error_message = f"Failed to start task: {str(e)}"
        db.commit()

        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error starting task: {str(e)}"
        )


@router.post("/{task_id}/cancel", response_model=MessageResponse)
async def cancel_task(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Cancel a running task

    Attempts to revoke the Celery task and marks it as cancelled
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    if task.status not in ["pending", "running"]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot cancel task with status: {task.status}"
        )

    # Revoke Celery task if exists
    if task.celery_task_id:
        try:
            from app.workers.celery_app import celery_app
            celery_app.control.revoke(task.celery_task_id, terminate=True)
        except Exception as e:
            print(f"Error revoking Celery task: {e}")

    # Update task status
    task.status = "cancelled"
    task.completed_at = datetime.utcnow()
    db.commit()

    return MessageResponse(message="Task cancelled successfully")


@router.get("/{task_id}/logs", response_model=TaskLogResponse)
async def get_task_logs(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Get task execution logs

    Returns the content of the task's log file if available
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    if not task.log_file_path:
        return TaskLogResponse(
            task_id=task.id,
            log_content="No logs available yet",
            log_file_path=None
        )

    # Read log file
    try:
        log_path = Path(task.log_file_path)
        if not log_path.exists():
            return TaskLogResponse(
                task_id=task.id,
                log_content="Log file not found",
                log_file_path=task.log_file_path
            )

        with open(log_path, "r") as f:
            log_content = f.read()

        return TaskLogResponse(
            task_id=task.id,
            log_content=log_content,
            log_file_path=task.log_file_path
        )

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error reading log file: {str(e)}"
        )


@router.delete("/{task_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_task(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Delete a task

    Only tasks with status 'pending', 'completed', 'failed', or 'cancelled' can be deleted.
    Running tasks must be cancelled first.
    """
    task = db.query(PipelineTask).join(Project).filter(
        PipelineTask.id == task_id,
        Project.user_id == current_user.id
    ).first()

    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Task not found"
        )

    if task.status == "running":
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot delete running task. Cancel it first."
        )

    # Delete log file if exists
    if task.log_file_path:
        try:
            log_path = Path(task.log_file_path)
            if log_path.exists():
                log_path.unlink()
        except Exception as e:
            print(f"Error deleting log file: {e}")

    # Delete task
    db.delete(task)
    db.commit()

    return None
