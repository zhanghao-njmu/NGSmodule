"""
System Job database model for tracking background operations
"""
from sqlalchemy import Column, String, Text, Float, DateTime, ForeignKey, JSON, Index
from sqlalchemy.orm import relationship
from app.core.types import UUID
import uuid

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class SystemJob(Base):
    """
    System job model for tracking background operations

    Tracks all system-level jobs including pipeline tasks,
    backups, cleanups, exports, and other administrative operations.
    """
    __tablename__ = "system_jobs"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)

    # Job classification
    type = Column(String(50), nullable=False, index=True)  # pipeline, backup, cleanup, export, import, system_task
    status = Column(String(20), nullable=False, default="pending", index=True)  # pending, running, completed, failed, cancelled

    # Job ownership
    user_id = Column(UUID(), ForeignKey("users.id"), nullable=True, index=True)
    username = Column(String(100), nullable=True)

    # Celery integration
    celery_task_id = Column(String(255), nullable=True, unique=True, index=True)

    # Progress tracking
    progress = Column(Float, default=0.0)  # 0.0 to 1.0
    message = Column(Text, nullable=True)

    # Job parameters and result
    parameters = Column(JSON, default=dict)
    result = Column(JSON, nullable=True)
    error = Column(Text, nullable=True)

    # Timestamps
    created_at = Column(DateTime, default=utc_now_naive, nullable=False, index=True)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Retry tracking
    retry_count = Column(String(10), default="0")
    parent_job_id = Column(UUID(), ForeignKey("system_jobs.id"), nullable=True)

    # Relationships
    user = relationship("User", foreign_keys=[user_id])
    parent_job = relationship("SystemJob", remote_side=[id])

    __table_args__ = (
        Index('ix_jobs_type_status', 'type', 'status'),
        Index('ix_jobs_user_created', 'user_id', 'created_at'),
    )

    def __repr__(self):
        return f"<SystemJob {self.type}/{self.status}: {self.id}>"

    def to_dict(self):
        """Convert job to dictionary"""
        return {
            "id": str(self.id),
            "type": self.type,
            "status": self.status,
            "user_id": str(self.user_id) if self.user_id else None,
            "username": self.username,
            "celery_task_id": self.celery_task_id,
            "progress": self.progress,
            "message": self.message,
            "parameters": self.parameters or {},
            "result": self.result,
            "error": self.error,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
            "retry_count": self.retry_count,
            "parent_job_id": str(self.parent_job_id) if self.parent_job_id else None,
        }
