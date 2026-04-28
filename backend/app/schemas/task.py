"""
Pipeline Task schemas
"""

from datetime import datetime
from typing import Any, Dict, List, Optional
from uuid import UUID

from pydantic import BaseModel, Field


class TaskBase(BaseModel):
    """Base task schema"""

    task_name: str = Field(..., min_length=1, max_length=100, description="Task name")
    task_type: Optional[str] = Field(None, description="Task type (e.g., 'RNA-seq', 'DNA-seq')")
    config: Dict[str, Any] = Field(default_factory=dict, description="Task configuration")


class TaskCreate(TaskBase):
    """Schema for creating a new task"""

    project_id: UUID = Field(..., description="Project ID to associate task with")


class TaskUpdate(BaseModel):
    """Schema for updating task"""

    task_name: Optional[str] = Field(None, min_length=1, max_length=100)
    task_type: Optional[str] = None
    status: Optional[str] = Field(None, pattern="^(pending|running|completed|failed|cancelled)$")
    progress: Optional[float] = Field(None, ge=0.0, le=100.0)
    error_message: Optional[str] = None
    config: Optional[Dict[str, Any]] = None


class TaskResponse(TaskBase):
    """Schema for task response"""

    id: UUID
    project_id: UUID
    status: str
    progress: float
    started_at: Optional[datetime]
    completed_at: Optional[datetime]
    error_message: Optional[str]
    celery_task_id: Optional[str]
    log_file_path: Optional[str]
    created_at: datetime

    model_config = {"from_attributes": True}


class TaskListResponse(BaseModel):
    """Schema for task list response"""

    total: int
    items: List[TaskResponse]


class TaskStats(BaseModel):
    """Schema for task statistics"""

    total_tasks: int
    pending_tasks: int
    running_tasks: int
    completed_tasks: int
    failed_tasks: int
    cancelled_tasks: int


class TaskLogResponse(BaseModel):
    """Schema for task log response"""

    task_id: UUID
    log_content: str
    log_file_path: Optional[str]


class TaskExecuteRequest(BaseModel):
    """Schema for executing a task"""

    pipeline_script: str = Field(..., description="Pipeline script name to execute")
    config: Dict[str, Any] = Field(default_factory=dict, description="Pipeline configuration parameters")
