"""
WebSocket endpoints for real-time updates
"""

import json
import logging
from datetime import datetime
from typing import Dict, Optional, Set

from fastapi import APIRouter, Query, WebSocket, WebSocketDisconnect

from app.core.database import SessionLocal
from app.core.deps import get_current_user_ws
from app.models.project import Project
from app.models.task import PipelineTask

logger = logging.getLogger(__name__)

router = APIRouter()


class ConnectionManager:
    """Manage WebSocket connections"""

    def __init__(self):
        # Store active connections: {user_id: {websocket, ...}}
        self.active_connections: Dict[str, Set[WebSocket]] = {}
        # Store task subscriptions: {task_id: {user_id, ...}}
        self.task_subscriptions: Dict[str, Set[str]] = {}

    async def connect(self, websocket: WebSocket, user_id: str):
        """Connect a new WebSocket client"""
        await websocket.accept()
        if user_id not in self.active_connections:
            self.active_connections[user_id] = set()
        self.active_connections[user_id].add(websocket)
        logger.info(f"WebSocket connected: user={user_id}")

    def disconnect(self, websocket: WebSocket, user_id: str):
        """Disconnect a WebSocket client"""
        if user_id in self.active_connections:
            self.active_connections[user_id].discard(websocket)
            if not self.active_connections[user_id]:
                del self.active_connections[user_id]
        logger.info(f"WebSocket disconnected: user={user_id}")

    def subscribe_to_task(self, task_id: str, user_id: str):
        """Subscribe user to task updates"""
        if task_id not in self.task_subscriptions:
            self.task_subscriptions[task_id] = set()
        self.task_subscriptions[task_id].add(user_id)

    def unsubscribe_from_task(self, task_id: str, user_id: str):
        """Unsubscribe user from task updates"""
        if task_id in self.task_subscriptions:
            self.task_subscriptions[task_id].discard(user_id)
            if not self.task_subscriptions[task_id]:
                del self.task_subscriptions[task_id]

    async def send_personal_message(self, message: dict, user_id: str):
        """Send message to specific user's connections"""
        if user_id in self.active_connections:
            disconnected = set()
            for websocket in self.active_connections[user_id]:
                try:
                    await websocket.send_json(message)
                except Exception as e:
                    logger.warning(f"Error sending to websocket: {e}")
                    disconnected.add(websocket)

            # Clean up disconnected websockets
            for ws in disconnected:
                self.active_connections[user_id].discard(ws)

    async def broadcast_task_update(self, task_id: str, message: dict):
        """Broadcast task update to all subscribed users"""
        if task_id in self.task_subscriptions:
            for user_id in self.task_subscriptions[task_id]:
                await self.send_personal_message(message, user_id)

    async def send_task_status(self, task_id: str, status: str, progress: float, message: str = ""):
        """Send task status update"""
        update_message = {
            "type": "task_update",
            "task_id": task_id,
            "status": status,
            "progress": progress,
            "message": message,
            "timestamp": datetime.utcnow().isoformat(),
        }
        await self.broadcast_task_update(task_id, update_message)


# Global connection manager
manager = ConnectionManager()


async def start_realtime_subscriber():
    """Hook called from the lifespan to start the Redis subscriber."""
    from app.services.realtime import RedisRealtimeSubscriber

    sub = RedisRealtimeSubscriber(manager)
    await sub.start()
    return sub


@router.websocket("/ws")
async def websocket_endpoint(websocket: WebSocket, token: Optional[str] = Query(None)):
    """
    WebSocket endpoint for real-time updates

    Connect with: ws://localhost:8000/api/v1/ws?token=<jwt_token>

    Message types from client:
    - {"type": "subscribe", "task_id": "uuid"} - Subscribe to task updates
    - {"type": "unsubscribe", "task_id": "uuid"} - Unsubscribe from task updates
    - {"type": "ping"} - Keep-alive ping

    Message types from server:
    - {"type": "task_update", "task_id": "uuid", "status": "running", "progress": 50.0, ...}
    - {"type": "pong"} - Response to ping
    - {"type": "error", "message": "error description"}
    """
    # Verify authentication
    if not token:
        await websocket.close(code=1008, reason="Missing authentication token")
        return

    try:
        # Verify user from token
        user = await get_current_user_ws(token)
        user_id = str(user.id)
    except Exception as e:
        await websocket.close(code=1008, reason=f"Authentication failed: {str(e)}")
        return

    # Connect
    await manager.connect(websocket, user_id)

    try:
        while True:
            # Receive message from client
            data = await websocket.receive_text()
            message = json.loads(data)

            message_type = message.get("type")

            if message_type == "subscribe":
                # Subscribe to task updates
                task_id = message.get("task_id")
                if task_id:
                    # Verify task belongs to user
                    db = SessionLocal()
                    try:
                        task = (
                            db.query(PipelineTask)
                            .join(Project)
                            .filter(PipelineTask.id == task_id, Project.user_id == user.id)
                            .first()
                        )

                        if task:
                            manager.subscribe_to_task(task_id, user_id)
                            await manager.send_personal_message(
                                {
                                    "type": "subscribed",
                                    "task_id": task_id,
                                    "status": task.status,
                                    "progress": task.progress,
                                },
                                user_id,
                            )
                        else:
                            await manager.send_personal_message(
                                {"type": "error", "message": "Task not found or access denied"}, user_id
                            )
                    finally:
                        db.close()

            elif message_type == "unsubscribe":
                # Unsubscribe from task updates
                task_id = message.get("task_id")
                if task_id:
                    manager.unsubscribe_from_task(task_id, user_id)
                    await manager.send_personal_message({"type": "unsubscribed", "task_id": task_id}, user_id)

            elif message_type == "ping":
                # Keep-alive ping
                await manager.send_personal_message(
                    {"type": "pong", "timestamp": datetime.utcnow().isoformat()}, user_id
                )

    except WebSocketDisconnect:
        manager.disconnect(websocket, user_id)
    except Exception as e:
        logger.error(f"WebSocket error: {e}")
        manager.disconnect(websocket, user_id)
