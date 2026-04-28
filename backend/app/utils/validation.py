"""
Validation utility functions
"""

from typing import Any, Optional, Type
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy.orm import Session


async def check_unique(
    db: Session,
    model: Type,
    field: str,
    value: Any,
    exclude_id: Optional[UUID] = None,
    error_msg: str = "Resource with this value already exists",
) -> None:
    """
    Check if a field value is unique in the database

    Args:
        db: Database session
        model: SQLAlchemy model class
        field: Field name to check
        value: Field value to check
        exclude_id: Record ID to exclude (useful for updates)
        error_msg: Custom error message

    Raises:
        HTTPException: If the value is not unique
    """
    # Build query to check existence
    query = db.query(model).filter(getattr(model, field) == value)

    # Exclude specific record (for updates)
    if exclude_id:
        query = query.filter(model.id != exclude_id)

    # Check if record exists
    existing = query.first()
    if existing:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)


async def check_unique_multi(
    db: Session,
    model: Type,
    fields: dict,
    exclude_id: Optional[UUID] = None,
    error_msg: str = "Resource with these values already exists",
) -> None:
    """
    Check if multiple field values are unique together

    Args:
        db: Database session
        model: SQLAlchemy model class
        fields: Dictionary of field names and values to check
        exclude_id: Record ID to exclude (useful for updates)
        error_msg: Custom error message

    Raises:
        HTTPException: If the combination is not unique
    """
    # Build query with all fields
    query = db.query(model)
    for field, value in fields.items():
        query = query.filter(getattr(model, field) == value)

    # Exclude specific record (for updates)
    if exclude_id:
        query = query.filter(model.id != exclude_id)

    # Check if record exists
    existing = query.first()
    if existing:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)


def validate_field_length(
    value: str, field_name: str, min_length: Optional[int] = None, max_length: Optional[int] = None
) -> None:
    """
    Validate string field length

    Args:
        value: String value to validate
        field_name: Name of the field (for error messages)
        min_length: Minimum length
        max_length: Maximum length

    Raises:
        HTTPException: If validation fails
    """
    if min_length and len(value) < min_length:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"{field_name} must be at least {min_length} characters"
        )

    if max_length and len(value) > max_length:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"{field_name} must be at most {max_length} characters"
        )


def validate_positive_number(value: int | float, field_name: str) -> None:
    """
    Validate that a number is positive

    Args:
        value: Number to validate
        field_name: Name of the field (for error messages)

    Raises:
        HTTPException: If validation fails
    """
    if value <= 0:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"{field_name} must be a positive number")


def validate_in_range(
    value: int | float,
    field_name: str,
    min_value: Optional[int | float] = None,
    max_value: Optional[int | float] = None,
) -> None:
    """
    Validate that a number is within a range

    Args:
        value: Number to validate
        field_name: Name of the field (for error messages)
        min_value: Minimum allowed value
        max_value: Maximum allowed value

    Raises:
        HTTPException: If validation fails
    """
    if min_value is not None and value < min_value:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"{field_name} must be at least {min_value}"
        )

    if max_value is not None and value > max_value:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"{field_name} must be at most {max_value}")
