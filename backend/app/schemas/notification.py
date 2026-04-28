"""
Notification schemas for API requests and responses
"""

from datetime import datetime
from typing import Any, Dict, Optional
from uuid import UUID

from pydantic import BaseModel, Field


class NotificationBase(BaseModel):
    """Base notification schema"""

    type: str = Field(..., description="Notification type", json_schema_extra={"example": "info"})
    title: str = Field(..., max_length=255, description="Notification title")
    message: str = Field(..., description="Notification message")
    data: Optional[Dict[str, Any]] = Field(None, description="Additional data")
    action_url: Optional[str] = Field(None, max_length=500, description="Action URL")
    priority: str = Field("normal", description="Priority level", pattern="^(low|normal|high|urgent)$")


class NotificationCreate(NotificationBase):
    """Schema for creating a notification"""

    user_id: UUID = Field(..., description="Target user ID")
    expires_at: Optional[datetime] = Field(None, description="Expiration datetime")


class NotificationUpdate(BaseModel):
    """Schema for updating a notification"""

    read: Optional[bool] = Field(None, description="Mark as read/unread")
    title: Optional[str] = Field(None, max_length=255)
    message: Optional[str] = None
    action_url: Optional[str] = Field(None, max_length=500)


class NotificationInDB(NotificationBase):
    """Schema for notification stored in database"""

    id: UUID
    user_id: UUID
    read: bool
    expires_at: Optional[datetime]
    created_at: datetime
    read_at: Optional[datetime]

    model_config = {"from_attributes": True}


class Notification(NotificationInDB):
    """Schema for notification API response"""


class NotificationList(BaseModel):
    """Schema for paginated notification list"""

    items: list[Notification]
    total: int
    page: int
    page_size: int
    unread_count: int


class NotificationSettingsBase(BaseModel):
    """Base notification settings schema"""

    email_enabled: bool = True
    email_task_completed: bool = True
    email_task_failed: bool = True
    email_system_alerts: bool = True
    app_enabled: bool = True
    app_task_updates: bool = True
    app_project_updates: bool = True
    app_system_alerts: bool = True
    push_enabled: bool = False


class NotificationSettingsCreate(NotificationSettingsBase):
    """Schema for creating notification settings"""

    user_id: UUID


class NotificationSettingsUpdate(BaseModel):
    """Schema for updating notification settings"""

    email_enabled: Optional[bool] = None
    email_task_completed: Optional[bool] = None
    email_task_failed: Optional[bool] = None
    email_system_alerts: Optional[bool] = None
    app_enabled: Optional[bool] = None
    app_task_updates: Optional[bool] = None
    app_project_updates: Optional[bool] = None
    app_system_alerts: Optional[bool] = None
    push_enabled: Optional[bool] = None


class NotificationSettingsInDB(NotificationSettingsBase):
    """Schema for notification settings in database"""

    id: UUID
    user_id: UUID
    updated_at: datetime

    model_config = {"from_attributes": True}


class NotificationSettings(NotificationSettingsInDB):
    """Schema for notification settings API response"""


class UnreadCount(BaseModel):
    """Schema for unread notification count"""

    count: int
    by_type: Dict[str, int] = Field(default_factory=dict)


class MarkAllReadResponse(BaseModel):
    """Schema for mark all as read response"""

    count: int = Field(..., description="Number of notifications marked as read")
    message: str = Field(..., description="Success message")


class WebSocketNotification(BaseModel):
    """Schema for WebSocket notification push"""

    event: str = Field("notification", description="Event type")
    data: Notification = Field(..., description="Notification data")
