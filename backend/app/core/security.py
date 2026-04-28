"""
Security utilities for authentication and authorization
"""

import base64
import hashlib
from datetime import timedelta
from functools import lru_cache
from typing import Any, Union

import bcrypt
import jwt
from cryptography.fernet import Fernet, InvalidToken

from app.core.config import settings
from app.core.datetime_utils import utc_now_naive

# bcrypt has a 72-byte limit. We use the bcrypt module directly to avoid
# passlib's compatibility issues with newer bcrypt versions (passlib's
# version detection breaks on bcrypt>=4 because `__about__` was removed).
_BCRYPT_MAX_BYTES = 72


def create_access_token(data: dict[str, Any], expires_delta: timedelta | None = None) -> str:
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
        expire = utc_now_naive() + timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)

    to_encode.update({"exp": expire})

    encoded_jwt = jwt.encode(to_encode, settings.JWT_SECRET_KEY, algorithm=settings.JWT_ALGORITHM)

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
        payload = jwt.decode(token, settings.JWT_SECRET_KEY, algorithms=[settings.JWT_ALGORITHM])
        return payload
    except jwt.PyJWTError:
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


# ----------- symmetric secret encryption (vendor credentials, API keys) -----------
#
# Stored secrets (e.g. lcbio email/password the user opted to save) are
# encrypted at rest with fernet. The key is derived from SECRET_KEY so
# rotating SECRET_KEY also invalidates stored ciphertexts — recommended
# behavior, since it forces re-entry of vendor passwords after a key
# rotation event.


@lru_cache(maxsize=1)
def _fernet() -> Fernet:
    # SHA-256 of SECRET_KEY → 32 bytes → base64-urlsafe → fernet key
    digest = hashlib.sha256(settings.SECRET_KEY.encode("utf-8")).digest()
    return Fernet(base64.urlsafe_b64encode(digest))


def encrypt_secret(plaintext: str) -> str:
    """Encrypt a string (e.g. vendor password) for at-rest storage."""
    return _fernet().encrypt(plaintext.encode("utf-8")).decode("utf-8")


def decrypt_secret(ciphertext: str) -> str:
    """Decrypt a fernet token. Raises InvalidToken on tamper / wrong key."""
    return _fernet().decrypt(ciphertext.encode("utf-8")).decode("utf-8")


def try_decrypt_secret(ciphertext: str) -> Union[str, None]:
    """Best-effort decrypt — returns None if the token is invalid / key
    has rotated. Useful when migrating or showing 'credentials need
    re-entry' messages."""
    try:
        return decrypt_secret(ciphertext)
    except (InvalidToken, ValueError):
        return None
