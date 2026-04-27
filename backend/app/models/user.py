"""
User model
"""
from sqlalchemy import Column, String, Boolean, DateTime, BigInteger
from app.core.types import UUID
from sqlalchemy.orm import relationship
import uuid
from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class User(Base):
    """User model for authentication and authorization"""

    __tablename__ = "users"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    username = Column(String(50), unique=True, nullable=False, index=True)
    email = Column(String(255), unique=True, nullable=False, index=True)
    password_hash = Column(String(255), nullable=False)
    full_name = Column(String(100))
    role = Column(String(20), default="user")  # 'user' or 'admin'
    organization = Column(String(100))
    is_active = Column(Boolean, default=True)
    storage_quota = Column(BigInteger, default=107374182400)  # 100GB in bytes
    storage_used = Column(BigInteger, default=0)
    last_login = Column(DateTime, nullable=True)
    created_at = Column(DateTime, default=utc_now_naive)
    updated_at = Column(DateTime, default=utc_now_naive, onupdate=utc_now_naive)

    # Relationships
    projects = relationship("Project", back_populates="owner", cascade="all, delete-orphan")
    notifications = relationship("Notification", back_populates="user", cascade="all, delete-orphan")
    notification_settings = relationship("NotificationSettings", back_populates="user", uselist=False, cascade="all, delete-orphan")

    def __repr__(self):
        return f"<User {self.username}>"

    @property
    def is_admin(self) -> bool:
        """Check if user is admin"""
        return self.role == "admin"

    @property
    def storage_available(self) -> int:
        """Get available storage in bytes"""
        return self.storage_quota - self.storage_used

    @property
    def storage_percent_used(self) -> float:
        """Get storage usage percentage"""
        if self.storage_quota == 0:
            return 0.0
        return (self.storage_used / self.storage_quota) * 100
