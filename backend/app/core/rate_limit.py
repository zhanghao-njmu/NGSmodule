"""
Rate limiting configuration for API endpoints

Uses slowapi for rate limiting with Redis storage
"""
import os
from slowapi import Limiter
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
from fastapi import Request, status
from fastapi.responses import JSONResponse
import logging

logger = logging.getLogger(__name__)

# Get Redis URL from environment for production, fallback to memory for development
# In production, set REDIS_URL environment variable (e.g., redis://localhost:6379/0)
RATE_LIMIT_STORAGE = os.getenv("REDIS_URL", "memory://")

# Initialize limiter with appropriate storage backend
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["100/minute"],  # Global default: 100 requests per minute
    storage_uri=RATE_LIMIT_STORAGE,
)

if RATE_LIMIT_STORAGE != "memory://":
    logger.info(f"Rate limiter using Redis storage")
else:
    logger.warning("Rate limiter using in-memory storage (not suitable for multi-instance deployment)")


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
    logger.warning(
        f"Rate limit exceeded for {get_remote_address(request)} "
        f"on {request.url.path}"
    )

    return JSONResponse(
        status_code=status.HTTP_429_TOO_MANY_REQUESTS,
        content={
            "error": "Rate limit exceeded",
            "detail": "Too many requests. Please slow down and try again later.",
            "code": "RATE_LIMIT_EXCEEDED",
            "retry_after": exc.detail  # Time until rate limit resets
        },
        headers={
            "Retry-After": str(exc.detail)
        }
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
