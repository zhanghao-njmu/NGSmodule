"""
Rate limiting configuration for API endpoints

Uses slowapi with Redis storage when available; falls back to in-memory
storage for local development/tests.
"""

import logging
import os

from fastapi import Request, status
from fastapi.responses import JSONResponse
from slowapi import Limiter
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

from app.core.config import settings

logger = logging.getLogger(__name__)


def _resolve_storage_uri() -> str:
    """Resolve rate limit storage from settings/env.

    Order of precedence:
      1. RATE_LIMIT_STORAGE env var (explicit override)
      2. settings.REDIS_URL (from .env / config) — only if Redis is reachable
      3. "memory://" (development fallback)

    During tests (pytest detected via sys.modules) we always use memory://.
    """
    import sys

    if "pytest" in sys.modules:
        return "memory://"

    explicit = os.getenv("RATE_LIMIT_STORAGE")
    if explicit:
        return explicit

    redis_url = getattr(settings, "REDIS_URL", None)
    if redis_url and redis_url not in ("", "memory://"):
        # Verify Redis is reachable; fall back to memory on connection error
        # so the app still starts when Redis is temporarily down.
        try:
            import redis as _redis

            _client = _redis.Redis.from_url(redis_url, socket_connect_timeout=1)
            _client.ping()
            _client.close()
            return redis_url
        except Exception as exc:
            logger.warning(f"Configured Redis ({redis_url}) is not reachable " f"({exc}); using in-memory rate limiter")

    return "memory://"


RATE_LIMIT_STORAGE = _resolve_storage_uri()

limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["100/minute"],
    storage_uri=RATE_LIMIT_STORAGE,
)

# Disable rate limiting under pytest. The default 100/min limit (and the
# tighter login/register limits) trip during test suites that hit dozens
# of endpoints in rapid succession with the same `testclient` IP.
import sys as _sys  # noqa: E402

if "pytest" in _sys.modules:
    limiter.enabled = False
    logger.info("Rate limiter disabled (pytest detected)")
elif RATE_LIMIT_STORAGE.startswith("redis"):
    logger.info(f"Rate limiter using Redis: {RATE_LIMIT_STORAGE}")
else:
    logger.warning("Rate limiter using in-memory storage " "(not suitable for multi-instance deployment)")


# Rate limit configurations for different endpoint types
class RateLimits:
    """
    Rate limit presets for different types of endpoints
    """

    # Authentication endpoints (more restrictive)
    LOGIN = "5/minute"  # 5 login attempts per minute
    REGISTER = "3/minute"  # 3 registration attempts per minute

    # File operations (moderate limits)
    FILE_UPLOAD = "10/minute"  # 10 file uploads per minute
    FILE_DOWNLOAD = "30/minute"  # 30 downloads per minute
    FILE_DELETE = "20/minute"  # 20 deletions per minute

    # CRUD operations (generous limits)
    CREATE = "30/minute"  # 30 creates per minute
    UPDATE = "50/minute"  # 50 updates per minute
    DELETE = "20/minute"  # 20 deletes per minute
    READ = "100/minute"  # 100 reads per minute (same as global)

    # Batch operations (more restrictive)
    BATCH_IMPORT = "5/minute"  # 5 batch imports per minute
    BATCH_EXPORT = "10/minute"  # 10 batch exports per minute

    # API-heavy operations
    STATS = "60/minute"  # 60 stats requests per minute
    SEARCH = "40/minute"  # 40 search requests per minute


def rate_limit_exceeded_handler(request: Request, exc: RateLimitExceeded) -> JSONResponse:
    """
    Custom error handler for rate limit exceeded

    Args:
        request: FastAPI request object
        exc: RateLimitExceeded exception

    Returns:
        JSONResponse with rate limit error details
    """
    logger.warning(f"Rate limit exceeded for {get_remote_address(request)} " f"on {request.url.path}")

    return JSONResponse(
        status_code=status.HTTP_429_TOO_MANY_REQUESTS,
        content={
            "error": "Rate limit exceeded",
            "detail": "Too many requests. Please slow down and try again later.",
            "code": "RATE_LIMIT_EXCEEDED",
            "retry_after": exc.detail,  # Time until rate limit resets
        },
        headers={"Retry-After": str(exc.detail)},
    )


def get_rate_limit_key(request: Request) -> str:
    """
    Generate rate limit key based on user or IP

    For authenticated requests, use user ID
    For anonymous requests, use IP address

    Args:
        request: FastAPI request object

    Returns:
        Rate limit key string
    """
    # Try to get user from request state (set by auth middleware)
    user = getattr(request.state, "user", None)

    if user and hasattr(user, "id"):
        return f"user:{user.id}"

    # Fallback to IP address
    return f"ip:{get_remote_address(request)}"


# Example: Configure Redis storage for production
def configure_redis_storage(redis_url: str):
    """
    Configure Redis as storage backend for rate limiting

    Args:
        redis_url: Redis connection URL (e.g., "redis://localhost:6379/0")

    Example:
        configure_redis_storage("redis://localhost:6379/1")
    """
    global limiter
    limiter = Limiter(
        key_func=get_remote_address,
        default_limits=["100/minute"],
        storage_uri=redis_url,
    )
    logger.info(f"Rate limiter configured with Redis storage: {redis_url}")
