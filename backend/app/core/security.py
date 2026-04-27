"""
Security utilities for authentication and authorization
"""
from datetime import timedelta
from typing import Any, Union
import bcrypt
from jose import jwt
from app.core.config import settings
from app.core.datetime_utils import utc_now_naive

# bcrypt has a 72-byte limit. We use the bcrypt module directly to avoid
# passlib's compatibility issues with newer bcrypt versions (passlib's
# version detection breaks on bcrypt>=4 because `__about__` was removed).
_BCRYPT_MAX_BYTES = 72


def create_access_token(
    data: dict[str, Any],
    expires_delta: timedelta | None = None
) -> str:
    """
    Create JWT access token

    Args:
        data: Payload data to encode
        expires_delta: Token expiration time

    Returns:
        Encoded JWT token
    """
    to_encode = data.copy()

    if expires_delta:
        expire = utc_now_naive() + expires_delta
    else:
        expire = utc_now_naive() + timedelta(
            minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES
        )

    to_encode.update({"exp": expire})

    encoded_jwt = jwt.encode(
        to_encode,
        settings.JWT_SECRET_KEY,
        algorithm=settings.JWT_ALGORITHM
    )

    return encoded_jwt


def verify_token(token: str) -> dict[str, Any] | None:
    """
    Verify and decode JWT token

    Args:
        token: JWT token to verify

    Returns:
        Decoded payload or None if invalid
    """
    try:
        payload = jwt.decode(
            token,
            settings.JWT_SECRET_KEY,
            algorithms=[settings.JWT_ALGORITHM]
        )
        return payload
    except jwt.JWTError:
        return None


def _truncate_password(password: str) -> bytes:
    """Encode and truncate password to bcrypt's 72-byte limit."""
    pw_bytes = password.encode("utf-8")
    if len(pw_bytes) > _BCRYPT_MAX_BYTES:
        pw_bytes = pw_bytes[:_BCRYPT_MAX_BYTES]
    return pw_bytes


def get_password_hash(password: str) -> str:
    """Hash a password with bcrypt."""
    salt = bcrypt.gensalt()
    hashed = bcrypt.hashpw(_truncate_password(password), salt)
    return hashed.decode("utf-8")


def verify_password(plain_password: str, hashed_password: str) -> bool:
    """Verify a password against its bcrypt hash."""
    try:
        return bcrypt.checkpw(
            _truncate_password(plain_password),
            hashed_password.encode("utf-8") if isinstance(hashed_password, str) else hashed_password,
        )
    except (ValueError, TypeError):
        return False
