"""
Real-time event broadcasting via Redis pub/sub.

Lets sync code (services, Celery tasks) publish events that get fanned out
to all FastAPI replicas, which in turn push them to connected WebSocket
clients.
"""

import asyncio
import json
import logging
from typing import Any, Dict, Optional

from app.core.config import settings

logger = logging.getLogger(__name__)


# Channel naming conventions
def user_channel(user_id: str) -> str:
    """Per-user notification channel."""
    return f"realtime:user:{user_id}"


def task_channel(task_id: str) -> str:
    """Per-task progress channel."""
    return f"realtime:task:{task_id}"


# ----------------------------------------------------------------------
# Publisher (sync, called from anywhere)
# ----------------------------------------------------------------------

_redis_publisher = None


def _get_publisher():
    """Lazy-initialize a sync Redis client used for publishing only."""
    global _redis_publisher
    if _redis_publisher is None:
        try:
            import redis
        except ImportError:
            logger.debug("redis package not installed; realtime events disabled")
            return None
        try:
            _redis_publisher = redis.Redis.from_url(settings.REDIS_URL, decode_responses=True)
            _redis_publisher.ping()
        except Exception as exc:
            logger.warning(f"Realtime publisher unavailable: {exc}")
            _redis_publisher = None
    return _redis_publisher


def publish_event(channel: str, event: Dict[str, Any]) -> bool:
    """Publish an event to a Redis channel. Returns True if delivered."""
    client = _get_publisher()
    if client is None:
        return False
    try:
        payload = json.dumps(event, default=str)
        client.publish(channel, payload)
        return True
    except Exception as exc:
        logger.warning(f"Failed to publish to {channel}: {exc}")
        return False


def publish_user_event(user_id: str, event_type: str, data: Dict[str, Any]) -> bool:
    """Publish an event to a specific user."""
    return publish_event(
        user_channel(str(user_id)),
        {"type": event_type, "data": data},
    )


def publish_task_event(task_id: str, event_type: str, data: Dict[str, Any]) -> bool:
    """Publish a task progress event."""
    return publish_event(
        task_channel(str(task_id)),
        {"type": event_type, "data": data},
    )


# ----------------------------------------------------------------------
# Subscriber (async, runs inside the FastAPI process)
# ----------------------------------------------------------------------


class RedisRealtimeSubscriber:
    """Subscribes to Redis pub/sub and forwards events to WebSocket clients.

    Maintains one Redis connection per process and uses pattern subscription
    so we can dynamically attach to new user/task channels without
    re-subscribing.
    """

    def __init__(self, manager):
        self.manager = manager  # WebSocket ConnectionManager
        self._task: Optional[asyncio.Task] = None
        self._stop = False

    async def start(self):
        if self._task is not None:
            return

        try:
            import redis.asyncio as redis_async
        except ImportError:
            logger.warning("redis package missing; realtime subscriber not starting")
            return

        self._stop = False
        self._task = asyncio.create_task(self._run(redis_async))
        logger.info("Realtime subscriber started")

    async def stop(self):
        self._stop = True
        if self._task:
            self._task.cancel()
            try:
                await self._task
            except asyncio.CancelledError:
                pass
            self._task = None

    async def _run(self, redis_async):
        try:
            client = redis_async.Redis.from_url(settings.REDIS_URL, decode_responses=True)
            pubsub = client.pubsub()
            await pubsub.psubscribe("realtime:user:*", "realtime:task:*")

            async for message in pubsub.listen():
                if self._stop:
                    break
                if message.get("type") not in ("pmessage", "message"):
                    continue
                channel = message.get("channel", "")
                try:
                    data = json.loads(message["data"])
                except (TypeError, ValueError):
                    continue

                if channel.startswith("realtime:user:"):
                    user_id = channel.split(":", 2)[2]
                    await self.manager.send_personal_message(data, user_id)
                elif channel.startswith("realtime:task:"):
                    task_id = channel.split(":", 2)[2]
                    await self.manager.broadcast_task_update(task_id, data)
        except asyncio.CancelledError:
            raise
        except Exception:
            logger.exception("Realtime subscriber crashed")
