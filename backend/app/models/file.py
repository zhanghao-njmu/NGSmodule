"""
File model
"""
from sqlalchemy import Column, String, DateTime, BigInteger, ForeignKey
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import relationship
import uuid
from datetime import datetime
from app.core.database import Base


class File(Base):
    """File model for uploaded files"""

    __tablename__ = "files"

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    sample_id = Column(UUID(as_uuid=True), ForeignKey("samples.id", ondelete="CASCADE"), nullable=False)
    filename = Column(String(255), nullable=False)
    file_path = Column(String(512), nullable=False)
    file_type = Column(String(20))  # 'fastq', 'bam', 'vcf', etc.
    file_size = Column(BigInteger)
    md5_checksum = Column(String(32))
    upload_status = Column(String(20), default="pending")  # 'pending', 'uploading', 'completed', 'failed'
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    sample = relationship("Sample", back_populates="files")

    def __repr__(self):
        return f"<File {self.filename}>"
