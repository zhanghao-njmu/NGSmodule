"""
Task Service - Business logic for pipeline task management
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy.orm import Session

from app.models.project import Project
from app.models.task import PipelineTask
from app.schemas.task import (
    TaskCreate,
    TaskExecuteRequest,
    TaskStats,
    TaskUpdate,
)
from app.workers.pipeline_tasks import run_ngs_pipeline

logger = logging.getLogger(__name__)


class TaskService:
    """Service class for task-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize TaskService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(self, task_id: UUID, user_id: UUID) -> Optional[PipelineTask]:
        """
        Get task by ID for a specific user

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Returns:
            Task if found and belongs to user, None otherwise
        """
        task = (
            self.db.query(PipelineTask)
            .join(Project)
            .filter(PipelineTask.id == task_id, Project.user_id == user_id)
            .first()
        )

        return task

    def get_by_id_or_raise(self, task_id: UUID, user_id: UUID) -> PipelineTask:
        """
        Get task by ID or raise 404 error

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Returns:
            Task

        Raises:
            HTTPException: If task not found or unauthorized
        """
        task = self.get_by_id(task_id, user_id)
        if not task:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Task {task_id} not found")
        return task

    def list_tasks(
        self,
        user_id: UUID,
        project_id: Optional[UUID] = None,
        status_filter: Optional[str] = None,
        task_type: Optional[str] = None,
    ) -> List[PipelineTask]:
        """
        List tasks for a user with optional filters

        Args:
            user_id: User UUID
            project_id: Optional project filter
            status_filter: Optional status filter
            task_type: Optional task type filter

        Returns:
            List of tasks
        """
        query = self.db.query(PipelineTask).join(Project).filter(Project.user_id == user_id)

        # Apply filters
        if project_id:
            query = query.filter(PipelineTask.project_id == project_id)

        if status_filter:
            query = query.filter(PipelineTask.status == status_filter)

        if task_type:
            query = query.filter(PipelineTask.task_type == task_type)

        tasks = query.order_by(PipelineTask.created_at.desc()).all()

        return tasks

    def get_stats(self, user_id: UUID, project_id: Optional[UUID] = None) -> TaskStats:
        """
        Get task statistics for a user

        Args:
            user_id: User UUID
            project_id: Optional project filter

        Returns:
            Task statistics
        """
        query = self.db.query(PipelineTask).join(Project).filter(Project.user_id == user_id)

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
            cancelled_tasks=cancelled_tasks,
        )

    # ============= CREATE OPERATIONS =============

    def create(self, user_id: UUID, task_data: TaskCreate) -> PipelineTask:
        """
        Create a new pipeline task

        Args:
            user_id: User UUID (for authorization)
            task_data: Task creation data

        Returns:
            Created task

        Raises:
            HTTPException: If project not found
        """
        # Verify project belongs to user
        project = self.db.query(Project).filter(Project.id == task_data.project_id, Project.user_id == user_id).first()

        if not project:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

        # Create task
        task = PipelineTask(
            project_id=task_data.project_id,
            task_name=task_data.task_name,
            task_type=task_data.task_type,
            config=task_data.config,
            status="pending",
        )

        self.db.add(task)
        self.db.commit()
        self.db.refresh(task)

        return task

    # ============= UPDATE OPERATIONS =============

    def update(self, task_id: UUID, user_id: UUID, update_data: TaskUpdate) -> PipelineTask:
        """
        Update task details

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)
            update_data: Update data

        Returns:
            Updated task

        Raises:
            HTTPException: If task not found
        """
        task = self.get_by_id_or_raise(task_id, user_id)

        # Update fields
        update_dict = update_data.model_dump(exclude_unset=True)
        for field, value in update_dict.items():
            setattr(task, field, value)

        self.db.commit()
        self.db.refresh(task)

        return task

    # ============= EXECUTION OPERATIONS =============

    def execute(self, task_id: UUID, user_id: UUID, execute_data: TaskExecuteRequest) -> Tuple[str, str]:
        """
        Execute a pipeline task

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)
            execute_data: Execution configuration

        Returns:
            Tuple of (message, celery_task_id)

        Raises:
            HTTPException: If task not found or already running/completed
        """
        task = self.get_by_id_or_raise(task_id, user_id)

        # Check if task is already running or completed
        if task.status in ["running", "completed"]:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Task is already {task.status}")

        # Update task configuration
        task.config = {**task.config, **execute_data.config}
        task.status = "pending"
        task.started_at = datetime.utcnow()
        self.db.commit()

        # Submit to Celery
        try:
            celery_task = run_ngs_pipeline.delay(
                task_id=str(task.id), pipeline_script=execute_data.pipeline_script, config=task.config
            )

            task.celery_task_id = celery_task.id
            self.db.commit()

            return (f"Task execution started with Celery ID: {celery_task.id}", celery_task.id)

        except Exception as e:
            task.status = "failed"
            task.error_message = f"Failed to start task: {str(e)}"
            self.db.commit()

            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Error starting task: {str(e)}"
            )

    def cancel(self, task_id: UUID, user_id: UUID) -> str:
        """
        Cancel a running task

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Returns:
            Success message

        Raises:
            HTTPException: If task not found or cannot be cancelled
        """
        task = self.get_by_id_or_raise(task_id, user_id)

        if task.status not in ["pending", "running"]:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail=f"Cannot cancel task with status: {task.status}"
            )

        # Revoke Celery task if exists
        if task.celery_task_id:
            try:
                from app.workers.celery_app import celery_app

                celery_app.control.revoke(task.celery_task_id, terminate=True)
            except Exception as e:
                logger.warning(f"Error revoking Celery task: {e}")

        # Update task status
        task.status = "cancelled"
        task.completed_at = datetime.utcnow()
        self.db.commit()

        return "Task cancelled successfully"

    # ============= LOG OPERATIONS =============

    def get_logs(self, task_id: UUID, user_id: UUID) -> Tuple[str, Optional[str]]:
        """
        Get task execution logs

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Returns:
            Tuple of (log_content, log_file_path)

        Raises:
            HTTPException: If task not found or error reading logs
        """
        task = self.get_by_id_or_raise(task_id, user_id)

        if not task.log_file_path:
            return ("No logs available yet", None)

        # Read log file
        try:
            log_path = Path(task.log_file_path)
            if not log_path.exists():
                return ("Log file not found", task.log_file_path)

            with open(log_path, "r") as f:
                log_content = f.read()

            return (log_content, task.log_file_path)

        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Error reading log file: {str(e)}"
            )

    # ============= DELETE OPERATIONS =============

    def delete(self, task_id: UUID, user_id: UUID) -> None:
        """
        Delete a task

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Raises:
            HTTPException: If task not found or is running
        """
        task = self.get_by_id_or_raise(task_id, user_id)

        if task.status == "running":
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail="Cannot delete running task. Cancel it first."
            )

        # Delete log file if exists
        if task.log_file_path:
            try:
                log_path = Path(task.log_file_path)
                if log_path.exists():
                    log_path.unlink()
            except Exception as e:
                logger.warning(f"Error deleting log file: {e}")

        # Delete task
        self.db.delete(task)
        self.db.commit()

    # ============= HELPER METHODS =============

    def get_tasks_by_project(
        self, project_id: UUID, user_id: UUID, status_filter: Optional[str] = None
    ) -> List[PipelineTask]:
        """
        Get all tasks for a specific project

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)
            status_filter: Optional status filter

        Returns:
            List of tasks

        Raises:
            HTTPException: If project not found
        """
        # Verify project belongs to user
        project = self.db.query(Project).filter(Project.id == project_id, Project.user_id == user_id).first()

        if not project:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

        query = self.db.query(PipelineTask).filter(PipelineTask.project_id == project_id)

        if status_filter:
            query = query.filter(PipelineTask.status == status_filter)

        tasks = query.order_by(PipelineTask.created_at.desc()).all()

        return tasks

    def update_progress(self, task_id: UUID, progress: float, status: Optional[str] = None) -> PipelineTask:
        """
        Update task progress (for use by workers)

        Args:
            task_id: Task UUID
            progress: Progress percentage (0-100)
            status: Optional status update

        Returns:
            Updated task

        Raises:
            HTTPException: If task not found
        """
        task = self.db.query(PipelineTask).filter(PipelineTask.id == task_id).first()

        if not task:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Task {task_id} not found")

        task.progress = progress

        if status:
            task.status = status

        if status == "completed":
            task.completed_at = datetime.utcnow()

        self.db.commit()
        self.db.refresh(task)

        return task
