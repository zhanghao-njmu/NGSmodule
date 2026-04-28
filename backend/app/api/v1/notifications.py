"""
Notifications API endpoints
"""

from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, Query, status
from sqlalchemy.orm import Session

from app.core.deps import get_current_active_user, get_db
from app.models.user import User
from app.schemas.notification import (
    MarkAllReadResponse,
    Notification,
    NotificationList,
)
from app.schemas.notification import NotificationSettings as NotificationSettingsSchema
from app.schemas.notification import (
    NotificationSettingsUpdate,
    UnreadCount,
)
from app.services.notification_service import NotificationService

router = APIRouter()


@router.get("", response_model=NotificationList)
async def get_notifications(
    skip: int = Query(0, ge=0, description="Number of records to skip"),
    limit: int = Query(20, ge=1, le=100, description="Maximum number of records"),
    unread_only: bool = Query(False, description="Show only unread notifications"),
    type: Optional[str] = Query(None, description="Filter by notification type"),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user),
):
    """
    Get user notifications with pagination

    - **skip**: Number of records to skip (for pagination)
    - **limit**: Maximum number of records to return
    - **unread_only**: Filter only unread notifications
    - **type**: Filter by notification type (info, warning, error, success)
    """
    service = NotificationService(db)

    notifications, total = service.get_notifications(
        user_id=current_user.id, skip=skip, limit=limit, unread_only=unread_only, notification_type=type
    )

    # Get unread count
    unread_count_obj = service.get_unread_count(current_user.id)

    return NotificationList(
        items=notifications, total=total, page=skip // limit + 1, page_size=limit, unread_count=unread_count_obj.count
    )


@router.get("/unread/count", response_model=UnreadCount)
async def get_unread_count(db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)):
    """
    Get count of unread notifications

    Returns total count and breakdown by notification type
    """
    service = NotificationService(db)
    return service.get_unread_count(current_user.id)


@router.get("/{notification_id}", response_model=Notification)
async def get_notification(
    notification_id: UUID, db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)
):
    """
    Get a specific notification by ID
    """
    service = NotificationService(db)
    notification = service.get_notification(notification_id, current_user.id)

    if not notification:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Notification not found")

    return notification


@router.put("/{notification_id}/read", response_model=Notification)
async def mark_notification_as_read(
    notification_id: UUID, db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)
):
    """
    Mark a notification as read
    """
    service = NotificationService(db)
    notification = service.mark_as_read(notification_id, current_user.id)

    if not notification:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Notification not found")

    return notification


@router.put("/read-all", response_model=MarkAllReadResponse)
async def mark_all_notifications_as_read(
    db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)
):
    """
    Mark all notifications as read
    """
    service = NotificationService(db)
    return service.mark_all_as_read(current_user.id)


@router.delete("/{notification_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_notification(
    notification_id: UUID, db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)
):
    """
    Delete a notification
    """
    service = NotificationService(db)
    deleted = service.delete_notification(notification_id, current_user.id)

    if not deleted:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Notification not found")

    return None


# ========== Notification Settings ==========


@router.get("/settings/current", response_model=NotificationSettingsSchema)
async def get_notification_settings(
    db: Session = Depends(get_db), current_user: User = Depends(get_current_active_user)
):
    """
    Get current user's notification settings
    """
    service = NotificationService(db)
    return service.get_settings(current_user.id)


@router.put("/settings/current", response_model=NotificationSettingsSchema)
async def update_notification_settings(
    settings_update: NotificationSettingsUpdate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_active_user),
):
    """
    Update notification settings

    Controls which types of notifications the user wants to receive
    """
    service = NotificationService(db)
    return service.update_settings(current_user.id, settings_update)
