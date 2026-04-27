"""
Pipeline Task model
"""
from sqlalchemy import Column, String, DateTime, Float, Text, ForeignKey
from app.core.types import UUID, JSONB
from sqlalchemy.orm import relationship
import uuid
from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class PipelineTask(Base):
    """Pipeline task model for tracking NGS analysis tasks"""

    __tablename__ = "pipeline_tasks"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    project_id = Column(UUID(), ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
    task_name = Column(String(100), nullable=False)
    task_type = Column(String(50))  # 'preAlignmentQC', 'Alignment', 'Quantification', etc.
    status = Column(String(20), default="pending")  # 'pending', 'running', 'completed', 'failed', 'cancelled'
    progress = Column(Float, default=0.0)  # 0-100
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    error_message = Column(Text)
    config = Column(JSONB, default={})  # Task configuration
    celery_task_id = Column(String(255))  # Celery task ID for tracking
    log_file_path = Column(String(512))
    created_at = Column(DateTime, default=utc_now_naive)

    # Relationships
    project = relationship("Project", back_populates="tasks")
    results = relationship("Result", back_populates="task", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<PipelineTask {self.task_name} - {self.status}>"
