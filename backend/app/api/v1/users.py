"""
User management API endpoints
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List
from app.core.database import get_db
from app.core.deps import get_current_user, get_current_admin
from app.models.user import User
from app.schemas.user import UserResponse, UserUpdate, UserAdminUpdate, UserStats, SystemStats
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.schemas.common import MessageResponse
from sqlalchemy import func

router = APIRouter()


@router.get("/me", response_model=UserResponse)
async def get_current_user_info(
    current_user: User = Depends(get_current_user)
):
    """
    Get current user information
    """
    return current_user


@router.put("/me", response_model=UserResponse)
async def update_current_user(
    user_update: UserUpdate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Update current user information

    - **full_name**: Optional full name
    - **organization**: Optional organization
    - **email**: Optional email
    """
    if user_update.full_name is not None:
        current_user.full_name = user_update.full_name

    if user_update.organization is not None:
        current_user.organization = user_update.organization

    if user_update.email is not None:
        # Check if email already exists
        existing_email = db.query(User).filter(
            User.email == user_update.email,
            User.id != current_user.id
        ).first()
        if existing_email:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Email already registered"
            )
        current_user.email = user_update.email

    db.commit()
    db.refresh(current_user)

    return current_user


@router.get("", response_model=List[UserResponse])
async def list_users(
    skip: int = 0,
    limit: int = 100,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    List all users (admin only)

    - **skip**: Number of users to skip
    - **limit**: Maximum number of users to return
    """
    users = db.query(User).offset(skip).limit(limit).all()
    return users


@router.get("/{user_id}", response_model=UserResponse)
async def get_user(
    user_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get user by ID (admin only)
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    return user


@router.put("/{user_id}", response_model=UserResponse)
async def update_user(
    user_id: str,
    user_update: UserAdminUpdate,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Update user by ID (admin only)

    Allows admin to update user role, quota, active status, etc.
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Update fields
    update_data = user_update.model_dump(exclude_unset=True)

    # Check email uniqueness if being updated
    if 'email' in update_data:
        existing_email = db.query(User).filter(
            User.email == update_data['email'],
            User.id != user_id
        ).first()
        if existing_email:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Email already registered"
            )

    for field, value in update_data.items():
        setattr(user, field, value)

    db.commit()
    db.refresh(user)

    return user


@router.post("/{user_id}/toggle", response_model=MessageResponse)
async def toggle_user_active_status(
    user_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Toggle user active status (admin only)
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Prevent admin from deactivating themselves
    if user.id == current_admin.id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot deactivate your own account"
        )

    user.is_active = not user.is_active
    db.commit()

    status_text = "activated" if user.is_active else "deactivated"
    return MessageResponse(message=f"User {user.username} {status_text} successfully")


@router.get("/{user_id}/stats", response_model=UserStats)
async def get_user_stats(
    user_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get user statistics (admin only)
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Count projects
    total_projects = db.query(func.count(Project.id)).filter(
        Project.user_id == user_id
    ).scalar() or 0

    # Count samples
    total_samples = db.query(func.count(Sample.id)).join(
        Project
    ).filter(Project.user_id == user_id).scalar() or 0

    # Count tasks
    total_tasks = db.query(func.count(PipelineTask.id)).join(
        Project
    ).filter(Project.user_id == user_id).scalar() or 0

    completed_tasks = db.query(func.count(PipelineTask.id)).join(
        Project
    ).filter(
        Project.user_id == user_id,
        PipelineTask.status == 'completed'
    ).scalar() or 0

    failed_tasks = db.query(func.count(PipelineTask.id)).join(
        Project
    ).filter(
        Project.user_id == user_id,
        PipelineTask.status == 'failed'
    ).scalar() or 0

    return UserStats(
        user_id=user.id,
        username=user.username,
        total_projects=total_projects,
        total_samples=total_samples,
        total_tasks=total_tasks,
        completed_tasks=completed_tasks,
        failed_tasks=failed_tasks,
        storage_used=user.storage_used,
        storage_quota=user.storage_quota,
        storage_percent=user.storage_percent_used
    )


@router.get("/stats/system", response_model=SystemStats)
async def get_system_stats(
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get system-wide statistics (admin only)
    """
    total_users = db.query(func.count(User.id)).scalar() or 0
    active_users = db.query(func.count(User.id)).filter(User.is_active == True).scalar() or 0

    total_projects = db.query(func.count(Project.id)).scalar() or 0
    total_samples = db.query(func.count(Sample.id)).scalar() or 0
    total_tasks = db.query(func.count(PipelineTask.id)).scalar() or 0

    running_tasks = db.query(func.count(PipelineTask.id)).filter(
        PipelineTask.status == 'running'
    ).scalar() or 0

    completed_tasks = db.query(func.count(PipelineTask.id)).filter(
        PipelineTask.status == 'completed'
    ).scalar() or 0

    failed_tasks = db.query(func.count(PipelineTask.id)).filter(
        PipelineTask.status == 'failed'
    ).scalar() or 0

    # Sum storage
    storage_stats = db.query(
        func.sum(User.storage_used).label('total_used'),
        func.sum(User.storage_quota).label('total_quota')
    ).first()

    return SystemStats(
        total_users=total_users,
        active_users=active_users,
        total_projects=total_projects,
        total_samples=total_samples,
        total_tasks=total_tasks,
        running_tasks=running_tasks,
        completed_tasks=completed_tasks,
        failed_tasks=failed_tasks,
        total_storage_used=storage_stats.total_used or 0,
        total_storage_quota=storage_stats.total_quota or 0
    )


@router.delete("/{user_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_user(
    user_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Delete user by ID (admin only)
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Prevent admin from deleting themselves
    if user.id == current_admin.id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot delete your own account"
        )

    db.delete(user)
    db.commit()

    return None
