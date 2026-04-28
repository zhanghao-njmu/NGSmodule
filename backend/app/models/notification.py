"""
Notification database model
"""

import uuid

from sqlalchemy import JSON, Boolean, Column, DateTime, ForeignKey, String, Text
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import UUID


class Notification(Base):
    """
    Notification model for user notifications

    Stores system notifications, task updates, and user alerts
    """

    __tablename__ = "notifications"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(), ForeignKey("users.id"), nullable=False, index=True)
    type = Column(String(50), nullable=False, index=True)  # info, warning, error, success, task_update
    title = Column(String(255), nullable=False)
    message = Column(Text, nullable=False)
    data = Column(JSON, nullable=True)  # Additional context data
    read = Column(Boolean, default=False, index=True)
    action_url = Column(String(500), nullable=True)  # Optional link to related resource
    priority = Column(String(20), default="normal")  # low, normal, high, urgent
    expires_at = Column(DateTime, nullable=True)  # Optional expiration date
    created_at = Column(DateTime, default=utc_now_naive, index=True)
    read_at = Column(DateTime, nullable=True)

    # Relationships
    user = relationship("User", back_populates="notifications")

    def __repr__(self):
        return f"<Notification {self.id}: {self.title}>"

    def to_dict(self):
        """Convert notification to dictionary"""
        return {
            "id": str(self.id),
            "user_id": str(self.user_id),
            "type": self.type,
            "title": self.title,
            "message": self.message,
            "data": self.data,
            "read": self.read,
            "action_url": self.action_url,
            "priority": self.priority,
            "expires_at": self.expires_at.isoformat() if self.expires_at else None,
            "created_at": self.created_at.isoformat(),
            "read_at": self.read_at.isoformat() if self.read_at else None,
        }


class NotificationSettings(Base):
    """
    User notification settings

    Controls which notifications a user wants to receive
    """

    __tablename__ = "notification_settings"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(), ForeignKey("users.id"), nullable=False, unique=True)

    # Email notifications
    email_enabled = Column(Boolean, default=True)
    email_task_completed = Column(Boolean, default=True)
    email_task_failed = Column(Boolean, default=True)
    email_system_alerts = Column(Boolean, default=True)

    # In-app notifications
    app_enabled = Column(Boolean, default=True)
    app_task_updates = Column(Boolean, default=True)
    app_project_updates = Column(Boolean, default=True)
    app_system_alerts = Column(Boolean, default=True)

    # Push notifications (future feature)
    push_enabled = Column(Boolean, default=False)

    updated_at = Column(DateTime, default=utc_now_naive, onupdate=utc_now_naive)

    # Relationships
    user = relationship("User", back_populates="notification_settings")

    def __repr__(self):
        return f"<NotificationSettings user_id={self.user_id}>"
