"""
Result model
"""
from sqlalchemy import Column, String, DateTime, Text, ForeignKey
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import relationship
import uuid
from datetime import datetime
from app.core.database import Base


class Result(Base):
    """Result model for storing analysis results"""

    __tablename__ = "results"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    task_id = Column(UUID(as_uuid=True), ForeignKey("pipeline_tasks.id", ondelete="CASCADE"), nullable=False)
    result_type = Column(String(50))  # 'qc_report', 'alignment', 'quantification', 'de_analysis'
    result_path = Column(Text)
    metadata = Column(JSONB, default={})  # Result metadata
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    task = relationship("PipelineTask", back_populates="results")

    def __repr__(self):
        return f"<Result {self.result_type}>"
