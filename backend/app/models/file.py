"""
File model
"""

import uuid

from sqlalchemy import BigInteger, Column, DateTime, ForeignKey, String
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import UUID


class File(Base):
    """File model for uploaded files"""

    __tablename__ = "files"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    sample_id = Column(UUID(), ForeignKey("samples.id", ondelete="CASCADE"), nullable=False)
    filename = Column(String(255), nullable=False)
    file_path = Column(String(512), nullable=False)
    file_type = Column(String(20))  # 'fastq', 'bam', 'vcf', etc.
    file_size = Column(BigInteger)
    md5_checksum = Column(String(32))
    upload_status = Column(String(20), default="pending")  # 'pending', 'uploading', 'completed', 'failed'
    created_at = Column(DateTime, default=utc_now_naive)

    # Relationships
    sample = relationship("Sample", back_populates="files")

    def __repr__(self):
        return f"<File {self.filename}>"
