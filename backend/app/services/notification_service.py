"""
Notification service for managing user notifications
"""

from datetime import datetime
from typing import List, Optional
from uuid import UUID

from sqlalchemy import and_, func
from sqlalchemy.orm import Session

from app.models.notification import Notification, NotificationSettings
from app.schemas.notification import (
    MarkAllReadResponse,
    NotificationCreate,
    NotificationSettingsUpdate,
    NotificationUpdate,
    UnreadCount,
)


class NotificationService:
    """Service for managing notifications"""

    def __init__(self, db: Session):
        self.db = db

    def create_notification(self, notification: NotificationCreate) -> Notification:
        """Create a new notification and broadcast it over WebSocket."""
        db_notification = Notification(**notification.model_dump())
        self.db.add(db_notification)
        self.db.commit()
        self.db.refresh(db_notification)

        # Publish to realtime channel so connected WebSocket clients get
        # immediate push updates (best-effort: failures are silent).
        try:
            from app.services.realtime import publish_user_event

            publish_user_event(
                str(db_notification.user_id),
                "notification",
                db_notification.to_dict(),
            )
        except Exception:
            pass

        return db_notification

    def get_notifications(
        self,
        user_id: UUID,
        skip: int = 0,
        limit: int = 20,
        unread_only: bool = False,
        notification_type: Optional[str] = None,
    ) -> tuple[List[Notification], int]:
        """
        Get user notifications with pagination

        Args:
            user_id: User ID
            skip: Number of records to skip
            limit: Maximum number of records to return
            unread_only: Filter only unread notifications
            notification_type: Filter by notification type

        Returns:
            Tuple of (notifications list, total count)
        """
        query = self.db.query(Notification).filter(Notification.user_id == user_id)

        if unread_only:
            query = query.filter(Notification.read == False)

        if notification_type:
            query = query.filter(Notification.type == notification_type)

        # Get total count
        total = query.count()

        # Get paginated results
        notifications = query.order_by(Notification.created_at.desc()).offset(skip).limit(limit).all()

        return notifications, total

    def get_notification(self, notification_id: UUID, user_id: UUID) -> Optional[Notification]:
        """
        Get a specific notification

        Args:
            notification_id: Notification ID
            user_id: User ID (for authorization)

        Returns:
            Notification or None
        """
        return (
            self.db.query(Notification)
            .filter(and_(Notification.id == notification_id, Notification.user_id == user_id))
            .first()
        )

    def update_notification(
        self, notification_id: UUID, user_id: UUID, notification_update: NotificationUpdate
    ) -> Optional[Notification]:
        """
        Update a notification

        Args:
            notification_id: Notification ID
            user_id: User ID (for authorization)
            notification_update: Update data

        Returns:
            Updated notification or None
        """
        db_notification = self.get_notification(notification_id, user_id)
        if not db_notification:
            return None

        update_data = notification_update.model_dump(exclude_unset=True)

        # If marking as read, set read_at timestamp
        if "read" in update_data and update_data["read"]:
            update_data["read_at"] = datetime.utcnow()
        elif "read" in update_data and not update_data["read"]:
            update_data["read_at"] = None

        for field, value in update_data.items():
            setattr(db_notification, field, value)

        self.db.commit()
        self.db.refresh(db_notification)
        return db_notification

    def mark_as_read(self, notification_id: UUID, user_id: UUID) -> Optional[Notification]:
        """
        Mark a notification as read

        Args:
            notification_id: Notification ID
            user_id: User ID

        Returns:
            Updated notification or None
        """
        update_data = NotificationUpdate(read=True)
        return self.update_notification(notification_id, user_id, update_data)

    def mark_all_as_read(self, user_id: UUID) -> MarkAllReadResponse:
        """
        Mark all user notifications as read

        Args:
            user_id: User ID

        Returns:
            Response with count of marked notifications
        """
        now = datetime.utcnow()
        result = (
            self.db.query(Notification)
            .filter(and_(Notification.user_id == user_id, Notification.read == False))
            .update({"read": True, "read_at": now}, synchronize_session=False)
        )
        self.db.commit()

        return MarkAllReadResponse(count=result, message=f"Marked {result} notifications as read")

    def delete_notification(self, notification_id: UUID, user_id: UUID) -> bool:
        """
        Delete a notification

        Args:
            notification_id: Notification ID
            user_id: User ID (for authorization)

        Returns:
            True if deleted, False if not found
        """
        db_notification = self.get_notification(notification_id, user_id)
        if not db_notification:
            return False

        self.db.delete(db_notification)
        self.db.commit()
        return True

    def get_unread_count(self, user_id: UUID) -> UnreadCount:
        """
        Get count of unread notifications

        Args:
            user_id: User ID

        Returns:
            UnreadCount with total and breakdown by type
        """
        total = (
            self.db.query(func.count(Notification.id))
            .filter(and_(Notification.user_id == user_id, Notification.read == False))
            .scalar()
            or 0
        )

        # Get count by type
        by_type_query = (
            self.db.query(Notification.type, func.count(Notification.id))
            .filter(and_(Notification.user_id == user_id, Notification.read == False))
            .group_by(Notification.type)
            .all()
        )

        by_type = {notification_type: count for notification_type, count in by_type_query}

        return UnreadCount(count=total, by_type=by_type)

    # ========== Notification Settings ==========

    def get_settings(self, user_id: UUID) -> NotificationSettings:
        """
        Get user notification settings

        Args:
            user_id: User ID

        Returns:
            NotificationSettings (creates default if not exists)
        """
        settings = self.db.query(NotificationSettings).filter(NotificationSettings.user_id == user_id).first()

        if not settings:
            # Create default settings
            settings = NotificationSettings(user_id=user_id)
            self.db.add(settings)
            self.db.commit()
            self.db.refresh(settings)

        return settings

    def update_settings(self, user_id: UUID, settings_update: NotificationSettingsUpdate) -> NotificationSettings:
        """
        Update notification settings

        Args:
            user_id: User ID
            settings_update: Settings update data

        Returns:
            Updated settings
        """
        settings = self.get_settings(user_id)

        update_data = settings_update.model_dump(exclude_unset=True)
        for field, value in update_data.items():
            setattr(settings, field, value)

        settings.updated_at = datetime.utcnow()
        self.db.commit()
        self.db.refresh(settings)
        return settings

    # ========== Notification Helpers ==========

    def create_task_notification(self, user_id: UUID, task_id: UUID, task_name: str, status: str) -> Notification:
        """
        Helper to create task status notification

        Args:
            user_id: User ID
            task_id: Task ID
            task_name: Task name
            status: Task status (completed, failed, etc.)

        Returns:
            Created notification
        """
        notification_types = {
            "completed": ("success", "Task Completed", f"Task '{task_name}' has completed successfully."),
            "failed": ("error", "Task Failed", f"Task '{task_name}' has failed. Please check the logs."),
            "running": ("info", "Task Started", f"Task '{task_name}' is now running."),
        }

        notif_type, title, message = notification_types.get(
            status, ("info", "Task Update", f"Task '{task_name}' status: {status}")
        )

        notification_data = NotificationCreate(
            user_id=user_id,
            type=notif_type,
            title=title,
            message=message,
            data={"task_id": str(task_id), "task_name": task_name, "status": status},
            action_url=f"/tasks/{task_id}",
            priority="normal" if status == "running" else "high",
        )

        return self.create_notification(notification_data)

    def create_system_notification(
        self, user_id: UUID, title: str, message: str, priority: str = "normal"
    ) -> Notification:
        """
        Helper to create system notification

        Args:
            user_id: User ID
            title: Notification title
            message: Notification message
            priority: Priority level

        Returns:
            Created notification
        """
        notification_data = NotificationCreate(
            user_id=user_id, type="system", title=title, message=message, priority=priority
        )

        return self.create_notification(notification_data)
