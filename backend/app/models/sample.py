"""
Sample model
"""
from sqlalchemy import Column, String, DateTime, ForeignKey
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import relationship
import uuid
from datetime import datetime
from app.core.database import Base


class Sample(Base):
    """Sample model for NGS samples"""

    __tablename__ = "samples"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    project_id = Column(UUID(as_uuid=True), ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
    sample_id = Column(String(100), nullable=False)
    run_id = Column(String(100))
    group_name = Column(String(50))
    layout = Column(String(10))  # 'PE' or 'SE'
    batch_id = Column(String(50))
    metadata = Column(JSONB, default={})  # Additional sample metadata
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    project = relationship("Project", back_populates="samples")
    files = relationship("File", back_populates="sample", cascade="all, delete-orphan")

    def __repr__(self):
        return f"<Sample {self.sample_id}>"
