"""
Admin API endpoints for user management, system configuration, and logs
"""
from fastapi import APIRouter, Depends, HTTPException, Query, status, Response
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from typing import Optional, List

from app.core.database import get_db
from app.core.security import get_current_admin_user
from app.models.user import User
from app.services.admin_service import AdminService
from app.schemas.admin import *


router = APIRouter()


# ============================================================================
# User Management Endpoints
# ============================================================================

@router.get("/users", response_model=UserListResponse)
async def get_all_users(
    skip: int = Query(0, ge=0),
    limit: int = Query(50, ge=1, le=100),
    role: Optional[UserRole] = None,
    is_active: Optional[bool] = None,
    search: Optional[str] = None,
    sort_by: str = Query("created_at", regex="^(created_at|username|email|role|storage_used)$"),
    sort_order: str = Query("desc", regex="^(asc|desc)$"),
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Get paginated list of all users (Admin only)

    **Parameters:**
    - skip: Number of records to skip (pagination)
    - limit: Maximum number of records to return
    - role: Filter by user role (user, admin)
    - is_active: Filter by active status
    - search: Search in username, email, full_name, organization
    - sort_by: Sort field (created_at, username, email, role, storage_used)
    - sort_order: Sort order (asc, desc)

    **Returns:**
    - Paginated list of users with basic info
    """
    service = AdminService(db)

    return service.get_users(
        skip=skip,
        limit=limit,
        role=role,
        is_active=is_active,
        search=search,
        sort_by=sort_by,
        sort_order=sort_order
    )


@router.get("/users/{user_id}", response_model=AdminUserDetail)
async def get_user_by_id(
    user_id: str,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Get detailed information about a specific user (Admin only)

    **Returns:**
    - Detailed user info including resource counts and activity
    """
    service = AdminService(db)

    user_detail = service.get_user_detail(user_id)
    if not user_detail:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    return user_detail


@router.put("/users/{user_id}", response_model=AdminUserDetail)
async def update_user_by_admin(
    user_id: str,
    update_data: UserUpdateRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Update user information (Admin only)

    **Parameters:**
    - username: New username (optional)
    - email: New email (optional)
    - full_name: New full name (optional)
    - organization: New organization (optional)
    - storage_quota: New storage quota in bytes (optional)

    **Returns:**
    - Updated user details
    """
    service = AdminService(db)

    # Check if user exists
    user = service.get_user_detail(user_id)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    # Update user
    updated_user = service.update_user(user_id, update_data)

    return updated_user


@router.put("/users/{user_id}/role", response_model=AdminUserDetail)
async def change_user_role(
    user_id: str,
    role_update: UserRoleUpdate,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Change user role (Admin only)

    **Parameters:**
    - role: New role (user, admin)

    **Returns:**
    - Updated user details

    **Security:**
    - Admins cannot demote themselves
    """
    service = AdminService(db)

    # Prevent admin from demoting themselves
    if str(current_admin.id) == user_id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="You cannot change your own role"
        )

    user = service.change_user_role(user_id, role_update.role)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    return user


@router.put("/users/{user_id}/activate", response_model=AdminUserDetail)
async def activate_deactivate_user(
    user_id: str,
    activation: UserActivationRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Activate or deactivate a user (Admin only)

    **Parameters:**
    - is_active: True to activate, False to deactivate
    - reason: Optional reason for the action

    **Returns:**
    - Updated user details

    **Security:**
    - Admins cannot deactivate themselves
    """
    service = AdminService(db)

    # Prevent admin from deactivating themselves
    if str(current_admin.id) == user_id and not activation.is_active:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="You cannot deactivate your own account"
        )

    user = service.activate_deactivate_user(user_id, activation.is_active, activation.reason)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    return user


@router.post("/users/{user_id}/reset-password", response_model=AdminOperationResponse)
async def reset_user_password(
    user_id: str,
    password_reset: PasswordResetRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Reset user password (Admin only)

    **Parameters:**
    - new_password: New password (min 8 chars, must contain uppercase, lowercase, digit)
    - notify_user: Send email notification to user (default: true)

    **Returns:**
    - Operation result

    **Security:**
    - Password is hashed before storage
    - Optional email notification to user
    """
    service = AdminService(db)

    success = service.reset_user_password(
        user_id,
        password_reset.new_password,
        password_reset.notify_user
    )

    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    return AdminOperationResponse(
        success=True,
        message="Password reset successfully",
        details={"notify_user": password_reset.notify_user}
    )


@router.delete("/users/{user_id}", response_model=AdminOperationResponse)
async def delete_user_by_admin(
    user_id: str,
    deletion: UserDeletionRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Delete a user account (Admin only)

    **Parameters:**
    - confirm: Must be true to proceed
    - transfer_data_to: Optional user ID to transfer data to
    - reason: Optional reason for deletion

    **Returns:**
    - Operation result

    **Security:**
    - Admins cannot delete themselves
    - Optionally transfer user's data to another user
    - Cascade deletion of related records

    **Warning:**
    - This action is irreversible!
    """
    service = AdminService(db)

    # Validate confirmation
    if not deletion.confirm:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Deletion must be confirmed"
        )

    # Prevent admin from deleting themselves
    if str(current_admin.id) == user_id:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="You cannot delete your own account"
        )

    success = service.delete_user(user_id, deletion.transfer_data_to)

    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )

    return AdminOperationResponse(
        success=True,
        message="User deleted successfully",
        details={
            "user_id": user_id,
            "data_transferred": deletion.transfer_data_to is not None,
            "reason": deletion.reason
        }
    )


# ============================================================================
# System Configuration Endpoints
# ============================================================================

@router.get("/config", response_model=SystemConfig)
async def get_system_configuration(
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Get complete system configuration (Admin only)

    **Returns:**
    - All system configuration categories:
      - general: App settings, quotas, session timeout
      - security: Auth settings, password policy
      - storage: Storage backend, retention policy
      - email: SMTP settings
      - notification: Notification preferences
      - pipeline: Pipeline execution settings
      - performance: API rate limits, caching
    """
    service = AdminService(db)

    return service.get_system_config()


@router.put("/config", response_model=SystemConfig)
async def update_system_configuration(
    config_update: ConfigUpdateRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Update system configuration (Admin only)

    **Parameters:**
    - category: Configuration category to update
    - updates: Dictionary of key-value pairs to update
    - reason: Optional reason for the change

    **Returns:**
    - Updated complete configuration

    **Note:**
    - Some changes may require server restart
    """
    service = AdminService(db)

    return service.update_system_config(
        category=config_update.category,
        updates=config_update.updates,
        admin_id=str(current_admin.id)
    )


@router.post("/config/reset", response_model=SystemConfig)
async def reset_system_configuration(
    reset_request: ConfigResetRequest,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Reset system configuration to defaults (Admin only)

    **Parameters:**
    - categories: List of categories to reset (None = reset all)
    - confirm: Must be true to proceed

    **Returns:**
    - Reset configuration

    **Warning:**
    - This will reset configuration to factory defaults!
    """
    if not reset_request.confirm:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Configuration reset must be confirmed"
        )

    service = AdminService(db)

    return service.reset_system_config(reset_request.categories)


# ============================================================================
# System Logs Endpoints
# ============================================================================

@router.get("/logs", response_model=LogResponse)
async def get_system_logs(
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    levels: Optional[List[LogLevel]] = Query(None),
    sources: Optional[List[LogSource]] = Query(None),
    search: Optional[str] = None,
    limit: int = Query(100, ge=1, le=1000),
    offset: int = Query(0, ge=0),
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Query system logs (Admin only)

    **Parameters:**
    - start_date: Filter logs from this date
    - end_date: Filter logs until this date
    - levels: Filter by log levels (debug, info, warning, error, critical)
    - sources: Filter by log sources (api, database, celery, pipeline, system, auth)
    - search: Search term in log messages
    - limit: Maximum number of logs to return (1-1000)
    - offset: Number of logs to skip (pagination)

    **Returns:**
    - List of log entries matching the criteria
    """
    service = AdminService(db)

    return service.get_logs(
        start_date=start_date,
        end_date=end_date,
        levels=levels,
        sources=sources,
        search=search,
        limit=limit,
        offset=offset
    )


@router.get("/logs/download")
async def download_system_logs(
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    levels: Optional[List[LogLevel]] = Query(None),
    sources: Optional[List[LogSource]] = Query(None),
    format: str = Query("json", regex="^(json|csv|txt)$"),
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Download system logs as a file (Admin only)

    **Parameters:**
    - start_date: Filter logs from this date
    - end_date: Filter logs until this date
    - levels: Filter by log levels
    - sources: Filter by log sources
    - format: File format (json, csv, txt)

    **Returns:**
    - Log file download
    """
    service = AdminService(db)

    log_file_path = service.download_logs(
        start_date=start_date,
        end_date=end_date,
        levels=levels,
        sources=sources,
        format=format
    )

    return FileResponse(
        path=log_file_path,
        filename=f"system_logs_{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.{format}",
        media_type="application/octet-stream"
    )


# ============================================================================
# System Health & Maintenance Endpoints
# ============================================================================

@router.get("/system/health", response_model=SystemHealth)
async def get_system_health(
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Get system health status (Admin only)

    **Returns:**
    - Overall system status
    - Individual service statuses:
      - PostgreSQL (database)
      - Redis (cache)
      - MinIO (object storage)
      - Celery (task queue)
    - Response times for each service
    - System uptime
    - System version

    **Status Values:**
    - healthy: All services operational
    - degraded: Some services have issues
    - down: Critical services are down
    """
    service = AdminService(db)

    return service.get_system_health()


@router.post("/system/cleanup", response_model=CleanupResponse)
async def cleanup_system(
    cleanup_options: CleanupOptions,
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Perform system cleanup (Admin only)

    **Parameters:**
    - old_logs: Delete logs older than retention period
    - temp_files: Delete temporary files
    - failed_tasks: Clean up failed task data
    - orphaned_files: Remove files not associated with any sample
    - old_notifications: Delete old read notifications
    - days_to_keep: Number of days to keep data (1-365)
    - dry_run: Preview what would be deleted without actually deleting

    **Returns:**
    - Cleanup results for each operation
    - Total items deleted
    - Total space freed (bytes)
    - Duration of operation

    **Recommendation:**
    - Run with dry_run=true first to preview changes
    """
    service = AdminService(db)

    return service.cleanup_system(cleanup_options)


# ============================================================================
# System Statistics (Admin)
# ============================================================================

@router.get("/system/stats", response_model=AdminSystemStats)
async def get_admin_system_statistics(
    current_admin: User = Depends(get_current_admin_user),
    db: Session = Depends(get_db)
):
    """
    Get comprehensive system statistics (Admin only)

    **Returns:**
    - User statistics (total, active, admin)
    - Resource counts (projects, samples, tasks, files)
    - Storage statistics (used, allocated, average per user)
    - Task statistics (today, this week, this month)
    - Success rates (today, this week, this month)
    - System metrics (CPU, memory, disk usage)
    - Active sessions

    **Use Case:**
    - System monitoring dashboard
    - Capacity planning
    - Performance analysis
    """
    service = AdminService(db)

    return service.get_admin_system_stats()


# ============================================================================
# Health Check
# ============================================================================

@router.get("/health")
async def admin_health_check(
    current_admin: User = Depends(get_current_admin_user)
):
    """
    Health check endpoint for admin service
    """
    return {
        "status": "healthy",
        "service": "admin",
        "timestamp": datetime.utcnow().isoformat()
    }
