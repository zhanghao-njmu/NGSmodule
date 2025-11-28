"""
Datetime utilities for timezone-aware datetime handling.

This module provides utilities for working with timezone-aware datetimes,
replacing the deprecated datetime.utcnow() pattern.
"""
from datetime import datetime, timezone


def utc_now() -> datetime:
    """
    Get current UTC time as a timezone-aware datetime.

    This replaces the deprecated datetime.utcnow() pattern.

    Returns:
        datetime: Current UTC time with timezone info
    """
    return datetime.now(timezone.utc)


def utc_now_naive() -> datetime:
    """
    Get current UTC time as a naive datetime (for SQLAlchemy compatibility).

    Some database drivers may not handle timezone-aware datetimes well.
    Use this when storing to databases that expect naive datetimes.

    Returns:
        datetime: Current UTC time without timezone info
    """
    return datetime.now(timezone.utc).replace(tzinfo=None)
