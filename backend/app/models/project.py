"""
Project model
"""
from sqlalchemy import Column, String, DateTime, Text, ForeignKey
from app.core.types import UUID, JSONB
from sqlalchemy.orm import relationship
import uuid
from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class Project(Base):
    """Project model for organizing NGS analysis"""

    __tablename__ = "projects"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(), ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    name = Column(String(100), nullable=False)
    description = Column(Text)
    project_type = Column(String(20))  # 'rna-seq', 'dna-seq', 'sc-rna-seq', etc.
    status = Column(String(20), default="active")  # 'active', 'archived', 'deleted'
    config = Column(JSONB, default={})  # Project configuration
    created_at = Column(DateTime, default=utc_now_naive)
    updated_at = Column(DateTime, default=utc_now_naive, onupdate=utc_now_naive)

    # Relationships
    owner = relationship("User", back_populates="projects")
    samples = relationship("Sample", back_populates="project", cascade="all, delete-orphan")
    tasks = relationship("PipelineTask", back_populates="project", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Project {self.name}>"
