"""
Data download model.

Represents a single delivery file fetched from a sequencing vendor
(联川 / 诺禾致源 / ...) into the local filesystem.
"""
import uuid
from sqlalchemy import Column, String, DateTime, BigInteger, Float, ForeignKey, Text
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import UUID


# Status state machine:
#   pending  -> queued, worker not yet picked up
#   running  -> lcbio Java service is downloading; progress_pct updated
#   completed -> file delivered fully; bytes_downloaded == file_size (if known)
#   failed   -> aborted with error_message
#   cancelled -> user-requested stop
DOWNLOAD_STATUS = ("pending", "running", "completed", "failed", "cancelled")


class DownloadJob(Base):
    """A single vendor → local download."""

    __tablename__ = "download_jobs"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(), ForeignKey("users.id", ondelete="CASCADE"), nullable=False)

    vendor = Column(String(32), nullable=False)  # "lc_bio" | "novogene" | ...
    source_path = Column(String(1024), nullable=False)  # vendor-side path (lcbio: obs_path)
    dest_path = Column(String(1024), nullable=False)  # local destination directory
    log_path = Column(String(1024))  # vendor download log (parsed for progress)

    status = Column(String(16), default="pending", nullable=False)
    progress_pct = Column(Float, default=0.0, nullable=False)
    bytes_downloaded = Column(BigInteger, default=0)
    file_size = Column(BigInteger)  # known after first progress event, optional
    error_message = Column(Text)

    started_at = Column(DateTime)
    finished_at = Column(DateTime)
    created_at = Column(DateTime, default=utc_now_naive, nullable=False)

    user = relationship("User")

    def __repr__(self):
        return f"<DownloadJob {self.vendor}:{self.source_path} status={self.status}>"
