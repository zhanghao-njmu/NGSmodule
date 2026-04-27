"""
Dependencies for dependency injection
"""
from uuid import UUID
from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from sqlalchemy.orm import Session
from jose import jwt, JWTError
from app.core.database import get_db
from app.core.config import settings
from app.core.security import verify_token
from app.core.permissions import verify_resource_ownership
from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.models.file import File as FileModel

oauth2_scheme = OAuth2PasswordBearer(tokenUrl=f"{settings.API_V1_PREFIX}/auth/login")


async def get_current_user(
    token: str = Depends(oauth2_scheme),
    db: Session = Depends(get_db)
) -> User:
    """
    Get current authenticated user

    Args:
        token: JWT access token
        db: Database session

    Returns:
        Current user

    Raises:
        HTTPException: If authentication fails
    """
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )

    payload = verify_token(token)
    if payload is None:
        raise credentials_exception

    user_id: str = payload.get("sub")
    if user_id is None:
        raise credentials_exception

    user = db.query(User).filter(User.id == user_id).first()
    if user is None:
        raise credentials_exception

    if not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User account is inactive"
        )

    return user


# Alias: existing API modules import this name, but get_current_user
# already enforces the is_active check. Keep both names for compatibility.
get_current_active_user = get_current_user


async def get_current_admin(
    current_user: User = Depends(get_current_user)
) -> User:
    """
    Get current admin user

    Args:
        current_user: Current authenticated user

    Returns:
        Current admin user

    Raises:
        HTTPException: If user is not admin
    """
    if not current_user.is_admin:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin access required"
        )

    return current_user


async def get_current_user_ws(token: str) -> User:
    """
    Get current authenticated user for WebSocket connections

    Args:
        token: JWT access token from query parameter

    Returns:
        Current user

    Raises:
        Exception: If authentication fails
    """
    from app.core.database import SessionLocal

    payload = verify_token(token)
    if payload is None:
        raise Exception("Invalid or expired token")

    user_id: str = payload.get("sub")
    if user_id is None:
        raise Exception("Invalid token payload")

    db = SessionLocal()
    try:
        user = db.query(User).filter(User.id == user_id).first()
        if user is None:
            raise Exception("User not found")

        if not user.is_active:
            raise Exception("User account is inactive")

        return user
    finally:
        db.close()


# Resource Ownership Verification Dependencies
# Refactored to use generic verify_resource_ownership function
async def get_user_project(
    project_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> Project:
    """
    Verify user ownership of a project

    Args:
        project_id: Project UUID
        current_user: Current authenticated user
        db: Database session

    Returns:
        Project if user has access

    Raises:
        HTTPException: If project not found or access denied
    """
    return verify_resource_ownership(
        db=db,
        model=Project,
        resource_id=project_id,
        current_user=current_user,
        resource_name="Project",
        filter_field="user_id"
    )


async def get_user_sample(
    sample_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> Sample:
    """
    Verify user ownership of a sample (through project)

    Args:
        sample_id: Sample UUID
        current_user: Current authenticated user
        db: Database session

    Returns:
        Sample if user has access

    Raises:
        HTTPException: If sample not found or access denied
    """
    return verify_resource_ownership(
        db=db,
        model=Sample,
        resource_id=sample_id,
        current_user=current_user,
        resource_name="Sample",
        join_models=[Project]
    )


async def get_user_task(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> PipelineTask:
    """
    Verify user ownership of a task (through project)

    Args:
        task_id: Task UUID
        current_user: Current authenticated user
        db: Database session

    Returns:
        Task if user has access

    Raises:
        HTTPException: If task not found or access denied
    """
    return verify_resource_ownership(
        db=db,
        model=PipelineTask,
        resource_id=task_id,
        current_user=current_user,
        resource_name="Task",
        join_models=[Project]
    )


async def get_user_file(
    file_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
) -> FileModel:
    """
    Verify user ownership of a file (through sample and project)

    Args:
        file_id: File UUID
        current_user: Current authenticated user
        db: Database session

    Returns:
        File if user has access

    Raises:
        HTTPException: If file not found or access denied
    """
    return verify_resource_ownership(
        db=db,
        model=FileModel,
        resource_id=file_id,
        current_user=current_user,
        resource_name="File",
        join_models=[Sample, Project]
    )
