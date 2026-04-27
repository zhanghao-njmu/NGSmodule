"""
Admin API endpoints for user management, system configuration, and logs
"""
from fastapi import APIRouter, Depends, HTTPException, Query, status, Response
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session
from typing import Optional, List

from app.core.database import get_db
from app.core.deps import get_current_admin
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
    current_admin: User = Depends(get_current_admin),
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
# Enhanced Features
# ============================================================================

@router.get("/system/metrics", response_model=SystemMetrics)
async def get_system_metrics(
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get detailed system metrics (Admin only)

    **Returns:**
    - CPU usage and load average
    - Memory usage (used, total, percentage)
    - Disk usage (used, total, percentage)
    - Network I/O statistics (optional)

    **Use Case:**
    - Real-time system monitoring
    - Performance analysis
    - Resource capacity planning
    """
    service = AdminService(db)
    return service.get_system_metrics()


@router.get("/alerts", response_model=AlertListResponse)
async def get_alerts(
    resolved: bool = False,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get system alerts (Admin only)

    **Parameters:**
    - resolved: Filter by resolution status (default: false, show unresolved)

    **Returns:**
    - List of alerts with type, severity, and details
    - Total count
    - Unresolved count

    **Alert Types:**
    - error: System errors requiring immediate attention
    - warning: Potential issues to monitor
    - info: Informational messages

    **Severity Levels:**
    - critical: Immediate action required
    - high: Important, address soon
    - medium: Monitor and plan
    - low: Informational
    """
    service = AdminService(db)
    return service.get_alerts(resolved=resolved)


@router.post("/alerts/{alert_id}/resolve")
async def resolve_alert(
    alert_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Resolve an alert (Admin only)

    **Parameters:**
    - alert_id: ID of the alert to resolve

    **Returns:**
    - Updated alert with resolved status
    - Resolution timestamp
    - Admin user who resolved it
    """
    service = AdminService(db)
    alert = service.resolve_alert(alert_id, str(current_admin.id))

    return {
        "success": True,
        "message": "Alert resolved successfully",
        "alert": alert
    }


@router.get("/audit-logs", response_model=List[AuditLogEntry])
async def get_audit_logs(
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    user_id: Optional[str] = None,
    action: Optional[str] = None,
    skip: int = 0,
    limit: int = 50,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get audit logs (Admin only)

    **Parameters:**
    - start_date: Filter by start date
    - end_date: Filter by end date
    - user_id: Filter by admin user ID
    - action: Filter by action type
    - skip: Pagination offset
    - limit: Number of records to return

    **Returns:**
    - List of audit log entries
    - Each entry includes:
      - Timestamp
      - Action performed
      - Admin user who performed it
      - Target user/resource
      - IP address
      - Details

    **Tracked Actions:**
    - User management (create, update, delete, role change)
    - System configuration changes
    - Password resets
    - Backups and restores
    - System cleanup operations
    """
    service = AdminService(db)
    return service.get_audit_logs(
        start_date=start_date,
        end_date=end_date,
        user_id=user_id,
        action=action,
        skip=skip,
        limit=limit
    )


@router.post("/audit-logs/export", response_model=ExportResult)
async def export_audit_logs(
    request: AuditLogExportRequest,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Export audit logs (Admin only)

    **Parameters:**
    - start_date: Start date for export range
    - end_date: End date for export range
    - user_id: Filter by admin user ID
    - action: Filter by action type
    - format: Export format (json or csv)

    **Returns:**
    - Download URL for the exported file
    - File name
    - File size
    - Expiration timestamp

    **Note:**
    - Export files expire after 24 hours
    - Large exports may take time to generate
    """
    service = AdminService(db)
    return service.export_audit_logs(
        start_date=request.start_date,
        end_date=request.end_date,
        user_id=request.user_id,
        action=request.action,
        format=request.format
    )


@router.get("/resources/usage", response_model=ResourceUsage)
async def get_resource_usage(
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Get resource usage summary (Admin only)

    **Returns:**
    - Storage: total and used disk space
    - Compute: active jobs and limits
    - Memory: total and used RAM

    **Use Case:**
    - Capacity planning
    - Resource allocation
    - Quota management
    """
    service = AdminService(db)
    return service.get_resource_usage()


@router.post("/backups", response_model=BackupInfo)
async def create_backup(
    request: BackupRequest,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Create system backup (Admin only)

    **Parameters:**
    - backup_type: Type of backup (full, incremental, database_only, files_only)
    - description: Optional description for the backup
    - compress: Whether to compress the backup (default: true)

    **Returns:**
    - Backup ID
    - File path
    - Size
    - Creation timestamp
    - Status

    **Backup Types:**
    - full: Complete system backup
    - incremental: Only changes since last backup
    - database_only: Database only
    - files_only: User files only

    **Note:**
    - Full backups are recommended weekly
    - Incremental backups can be done daily
    """
    service = AdminService(db)
    return service.create_backup(
        backup_type=request.backup_type,
        description=request.description,
        admin_user_id=str(current_admin.id),
        compress=request.compress
    )


@router.get("/backups", response_model=BackupListResponse)
async def list_backups(
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    List all backups (Admin only)

    **Returns:**
    - List of all system backups
    - Each backup includes:
      - ID
      - Type
      - File path
      - Size
      - Creation date
      - Status

    **Use Case:**
    - Backup management
    - Restore planning
    - Storage cleanup
    """
    service = AdminService(db)
    return service.list_backups()


@router.get("/jobs", response_model=JobListResponse)
async def list_jobs(
    status: Optional[JobStatus] = None,
    skip: int = 0,
    limit: int = 50,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    List system jobs (Admin only)

    **Parameters:**
    - status: Filter by job status (pending, running, completed, failed, cancelled)
    - skip: Pagination offset
    - limit: Number of records to return

    **Returns:**
    - List of jobs with details
    - Each job includes:
      - ID
      - Type (pipeline, backup, cleanup, export, import)
      - Status
      - User
      - Timestamps
      - Progress
      - Result or error message

    **Use Case:**
    - Job monitoring
    - Performance analysis
    - Error diagnosis
    """
    service = AdminService(db)
    return service.list_jobs(status=status, skip=skip, limit=limit)


@router.post("/jobs/{job_id}/cancel")
async def cancel_job(
    job_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Cancel a running job (Admin only)

    **Parameters:**
    - job_id: ID of the job to cancel

    **Returns:**
    - Updated job with cancelled status
    - Cancellation message

    **Note:**
    - Only pending or running jobs can be cancelled
    - Cancellation may take time to complete
    """
    service = AdminService(db)
    job = service.cancel_job(job_id)

    return {
        "success": True,
        "message": "Job cancelled successfully",
        "job": job
    }


@router.post("/jobs/{job_id}/retry")
async def retry_job(
    job_id: str,
    current_admin: User = Depends(get_current_admin),
    db: Session = Depends(get_db)
):
    """
    Retry a failed job (Admin only)

    **Parameters:**
    - job_id: ID of the job to retry

    **Returns:**
    - New job with pending status
    - Retry message

    **Note:**
    - Only failed jobs can be retried
    - A new job will be created with the same parameters
    """
    service = AdminService(db)
    job = service.retry_job(job_id)

    return {
        "success": True,
        "message": "Job queued for retry",
        "job": job
    }


# ============================================================================
# Health Check
# ============================================================================

@router.get("/health")
async def admin_health_check(
    current_admin: User = Depends(get_current_admin)
):
    """
    Health check endpoint for admin service
    """
    return {
        "status": "healthy",
        "service": "admin",
        "timestamp": datetime.utcnow().isoformat()
    }
