"""
System Backup database model
"""

import uuid

from sqlalchemy import JSON, BigInteger, Boolean, Column, DateTime, ForeignKey, Index, String, Text
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import UUID


class SystemBackup(Base):
    """
    System backup model

    Tracks all system backups including database dumps, file backups,
    and full system snapshots.
    """

    __tablename__ = "system_backups"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)

    # Backup classification
    backup_type = Column(String(50), nullable=False, index=True)  # full, incremental, database_only, files_only
    status = Column(
        String(20), nullable=False, default="pending", index=True
    )  # pending, in_progress, completed, failed

    # Backup file info
    file_path = Column(String(1024), nullable=True)
    file_name = Column(String(255), nullable=True)
    size = Column(BigInteger, default=0)  # bytes
    compressed = Column(Boolean, default=True)
    checksum = Column(String(128), nullable=True)  # SHA256

    # Backup details
    description = Column(Text, nullable=True)
    backup_metadata = Column(JSON, default=dict)
    error_message = Column(Text, nullable=True)

    # Ownership
    created_by = Column(UUID(), ForeignKey("users.id"), nullable=False)

    # Timestamps
    created_at = Column(DateTime, default=utc_now_naive, nullable=False, index=True)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    expires_at = Column(DateTime, nullable=True, index=True)  # Retention policy

    # Relationships
    creator = relationship("User", foreign_keys=[created_by])

    __table_args__ = (
        Index("ix_backups_type_status", "backup_type", "status"),
        Index("ix_backups_created_by_at", "created_by", "created_at"),
    )

    def __repr__(self):
        return f"<SystemBackup {self.backup_type}/{self.status}: {self.file_name}>"

    def to_dict(self):
        """Convert backup to dictionary"""
        return {
            "id": str(self.id),
            "backup_type": self.backup_type,
            "status": self.status,
            "file_path": self.file_path,
            "file_name": self.file_name,
            "size": self.size,
            "compressed": self.compressed,
            "checksum": self.checksum,
            "description": self.description,
            "metadata": self.backup_metadata or {},
            "error_message": self.error_message,
            "created_by": str(self.created_by),
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
            "expires_at": self.expires_at.isoformat() if self.expires_at else None,
        }
