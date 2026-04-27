"""
Sample model
"""
from sqlalchemy import Column, String, DateTime, ForeignKey
from app.core.types import UUID, JSONB
from sqlalchemy.orm import relationship
import uuid
from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class Sample(Base):
    """Sample model for NGS samples"""

    __tablename__ = "samples"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    project_id = Column(UUID(), ForeignKey("projects.id", ondelete="CASCADE"), nullable=False)
    sample_id = Column(String(100), nullable=False)
    run_id = Column(String(100))
    group_name = Column(String(50))
    layout = Column(String(10))  # 'PE' or 'SE'
    batch_id = Column(String(50))
    # NOTE: Python attribute renamed to avoid clash with SQLAlchemy's reserved `metadata` name.
    # Database column is still named `metadata` for backward compatibility.
    sample_metadata = Column("metadata", JSONB, default={})
    created_at = Column(DateTime, default=utc_now_naive)

    # Relationships
    project = relationship("Project", back_populates="samples")
    files = relationship("File", back_populates="sample", cascade="all, delete-orphan")

    @property
    def metadata_dict(self) -> dict:
        """Backward-compatible accessor for sample metadata."""
        return self.sample_metadata or {}

    def __repr__(self):
        return f"<Sample {self.sample_id}>"
