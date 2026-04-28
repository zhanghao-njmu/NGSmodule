"""
Generic permission and resource ownership verification utilities

Provides reusable functions for checking user access to resources
"""

from typing import Optional, Type, TypeVar
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy.orm import Query, Session

from app.models.project import Project
from app.models.user import User

# Generic type for database models
ModelType = TypeVar("ModelType")


def verify_resource_ownership(
    db: Session,
    model: Type[ModelType],
    resource_id: UUID,
    current_user: User,
    resource_name: str = "Resource",
    join_models: Optional[list] = None,
    filter_field: Optional[str] = None,
) -> ModelType:
    """
    Generic function to verify user ownership of a resource

    Args:
        db: Database session
        model: SQLAlchemy model class
        resource_id: UUID of the resource
        current_user: Current authenticated user
        resource_name: Name of resource for error messages (e.g., "Project", "Sample")
        join_models: List of models to join (for nested ownership checks)
        filter_field: Field name to filter by user (defaults to "user_id" via Project)

    Returns:
        The resource if user has access

    Raises:
        HTTPException: 404 if resource not found or access denied

    Examples:
        # Direct ownership (user_id field):
        project = verify_resource_ownership(
            db, Project, project_id, current_user, "Project"
        )

        # Nested ownership (through Project):
        sample = verify_resource_ownership(
            db, Sample, sample_id, current_user, "Sample",
            join_models=[Project]
        )
    """
    # Admin can access all resources
    if current_user.is_admin:
        resource = db.query(model).filter(model.id == resource_id).first()
    else:
        # Build query with joins if needed
        query = db.query(model)

        if join_models:
            for join_model in join_models:
                query = query.join(join_model)

        # Apply filter
        if join_models and Project in join_models:
            # Filter by project user_id for nested resources
            query = query.filter(model.id == resource_id, Project.user_id == current_user.id)
        else:
            # Direct ownership - assume model has user_id
            query = query.filter(model.id == resource_id, getattr(model, filter_field or "user_id") == current_user.id)

        resource = query.first()

    if not resource:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"{resource_name} not found or access denied")

    return resource


def require_role(role: str):
    """
    Decorator that requires the current user to have a specific role.

    Used as `@require_role("admin")` on FastAPI endpoints. Expects the
    decorated coroutine to receive `current_user: User` as a kwarg.
    """
    from functools import wraps

    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            current_user = kwargs.get("current_user")
            if current_user is None:
                # Locate User instance in positional args as fallback
                for arg in args:
                    if isinstance(arg, User):
                        current_user = arg
                        break

            if current_user is None:
                raise HTTPException(
                    status_code=status.HTTP_401_UNAUTHORIZED,
                    detail="Authentication required",
                )

            user_role = getattr(current_user, "role", None)
            if user_role != role and not (role == "admin" and getattr(current_user, "is_admin", False)):
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail=f"{role.title()} access required",
                )

            return await func(*args, **kwargs)

        return wrapper

    return decorator


def require_admin(current_user: User) -> None:
    """
    Verify that the current user has admin privileges

    Args:
        current_user: Current authenticated user

    Raises:
        HTTPException: 403 if user is not admin
    """
    if not current_user.is_admin:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Admin access required")


def require_active_user(current_user: User) -> None:
    """
    Verify that the current user account is active

    Args:
        current_user: Current authenticated user

    Raises:
        HTTPException: 403 if user account is inactive
    """
    if not current_user.is_active:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="User account is inactive")


def check_resource_quota(
    current_user: User, quota_field: str, limit_field: str, resource_name: str = "resource"
) -> None:
    """
    Check if user has quota available for creating resources

    Args:
        current_user: Current authenticated user
        quota_field: Field name for current usage (e.g., "project_count")
        limit_field: Field name for limit (e.g., "max_projects")
        resource_name: Name of resource for error message

    Raises:
        HTTPException: 403 if quota exceeded

    Example:
        check_resource_quota(
            current_user,
            "project_count",
            "max_projects",
            "project"
        )
    """
    current = getattr(current_user, quota_field, 0)
    limit = getattr(current_user, limit_field, float("inf"))

    if current >= limit:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"{resource_name.capitalize()} quota exceeded. " f"Current: {current}, Limit: {limit}",
        )


def build_user_filter(
    query: Query, current_user: User, ownership_model: Type = None, ownership_field: str = "user_id"
) -> Query:
    """
    Add user ownership filter to a query (admin sees all)

    Args:
        query: SQLAlchemy query object
        current_user: Current authenticated user
        ownership_model: Model containing the ownership field (if different from queried model)
        ownership_field: Field name for ownership (default: "user_id")

    Returns:
        Filtered query

    Example:
        # Direct ownership:
        query = db.query(Project)
        query = build_user_filter(query, current_user)

        # Nested ownership:
        query = db.query(Sample).join(Project)
        query = build_user_filter(query, current_user, Project)
    """
    if current_user.is_admin:
        return query

    model = ownership_model if ownership_model else query.column_descriptions[0]["type"]
    return query.filter(getattr(model, ownership_field) == current_user.id)
