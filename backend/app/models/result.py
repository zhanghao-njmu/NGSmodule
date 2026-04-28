"""
Result model
"""

import uuid

from sqlalchemy import Column, DateTime, ForeignKey, String, Text
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import JSONB, UUID


class Result(Base):
    """Result model for storing analysis results"""

    __tablename__ = "results"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    task_id = Column(UUID(), ForeignKey("pipeline_tasks.id", ondelete="CASCADE"), nullable=False)
    result_type = Column(String(50))  # 'qc_report', 'alignment', 'quantification', 'de_analysis'
    result_path = Column(Text)
    # NOTE: Python attribute renamed to avoid clash with SQLAlchemy's reserved `metadata` name.
    # Database column is still named `metadata` for backward compatibility.
    result_metadata = Column("metadata", JSONB, default={})
    created_at = Column(DateTime, default=utc_now_naive)

    # Relationships
    task = relationship("PipelineTask", back_populates="results")

    def __repr__(self):
        return f"<Result {self.result_type}>"
