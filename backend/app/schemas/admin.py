"""
Admin schemas for user management, system configuration, and logs
"""
from pydantic import BaseModel, EmailStr, Field, validator
from typing import List, Dict, Optional, Any
from datetime import datetime
from enum import Enum


# ============================================================================
# User Management
# ============================================================================

class UserRole(str, Enum):
    """User role enumeration"""
    USER = "user"
    ADMIN = "admin"


class UserStatus(str, Enum):
    """User status enumeration"""
    ACTIVE = "active"
    INACTIVE = "inactive"
    SUSPENDED = "suspended"


class AdminUserList(BaseModel):
    """User in admin list view"""
    id: str
    username: str
    email: EmailStr
    full_name: Optional[str] = None
    role: UserRole
    organization: Optional[str] = None
    is_active: bool
    storage_used: int
    storage_quota: int
    storage_percent: float
    last_login: Optional[datetime] = None
    created_at: datetime
    updated_at: datetime

    class Config:
        from_attributes = True


class AdminUserDetail(AdminUserList):
    """Detailed user information for admin"""
    total_projects: int = 0
    total_samples: int = 0
    total_tasks: int = 0
    total_files: int = 0
    last_activity: Optional[datetime] = None
    login_count: Optional[int] = 0


class UserListResponse(BaseModel):
    """Paginated user list response"""
    users: List[AdminUserList]
    total: int
    page: int
    page_size: int
    total_pages: int


class UserUpdateRequest(BaseModel):
    """Admin update user request"""
    username: Optional[str] = Field(None, min_length=3, max_length=50)
    email: Optional[EmailStr] = None
    full_name: Optional[str] = Field(None, max_length=100)
    organization: Optional[str] = Field(None, max_length=100)
    storage_quota: Optional[int] = Field(None, gt=0)

    @validator('username')
    def validate_username(cls, v):
        if v and not v.replace('_', '').replace('-', '').isalnum():
            raise ValueError('Username must contain only letters, numbers, underscores, and hyphens')
        return v


class UserRoleUpdate(BaseModel):
    """Update user role"""
    role: UserRole


class UserActivationRequest(BaseModel):
    """Activate or deactivate user"""
    is_active: bool
    reason: Optional[str] = None


class PasswordResetRequest(BaseModel):
    """Admin reset user password"""
    new_password: str = Field(..., min_length=8, max_length=100)
    notify_user: bool = True  # Send email notification to user

    @validator('new_password')
    def validate_password(cls, v):
        if not any(c.isupper() for c in v):
            raise ValueError('Password must contain at least one uppercase letter')
        if not any(c.islower() for c in v):
            raise ValueError('Password must contain at least one lowercase letter')
        if not any(c.isdigit() for c in v):
            raise ValueError('Password must contain at least one digit')
        return v


class UserDeletionRequest(BaseModel):
    """Delete user request"""
    confirm: bool = True
    transfer_data_to: Optional[str] = None  # Transfer user's data to another user
    reason: Optional[str] = None


class BulkUserOperation(BaseModel):
    """Bulk operation on multiple users"""
    user_ids: List[str] = Field(..., min_length=1)
    operation: str  # activate, deactivate, delete, change_role
    parameters: Optional[Dict[str, Any]] = None


# ============================================================================
# System Configuration
# ============================================================================

class SystemConfigCategory(str, Enum):
    """System configuration category"""
    GENERAL = "general"
    SECURITY = "security"
    STORAGE = "storage"
    EMAIL = "email"
    NOTIFICATION = "notification"
    PIPELINE = "pipeline"
    PERFORMANCE = "performance"


class ConfigItem(BaseModel):
    """Single configuration item"""
    key: str
    value: Any
    category: SystemConfigCategory
    description: Optional[str] = None
    value_type: str  # string, int, bool, float, json
    default_value: Any
    is_secret: bool = False  # Hide value in responses
    requires_restart: bool = False
    last_modified: Optional[datetime] = None
    modified_by: Optional[str] = None


class SystemConfig(BaseModel):
    """Complete system configuration"""
    general: Dict[str, Any]
    security: Dict[str, Any]
    storage: Dict[str, Any]
    email: Dict[str, Any]
    notification: Dict[str, Any]
    pipeline: Dict[str, Any]
    performance: Dict[str, Any]
    last_updated: datetime


class ConfigUpdateRequest(BaseModel):
    """Update system configuration"""
    category: SystemConfigCategory
    updates: Dict[str, Any]
    reason: Optional[str] = None


class ConfigResetRequest(BaseModel):
    """Reset configuration to defaults"""
    categories: Optional[List[SystemConfigCategory]] = None  # None = reset all
    confirm: bool = True


# ============================================================================
# System Logs
# ============================================================================

class LogLevel(str, Enum):
    """Log level enumeration"""
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class LogSource(str, Enum):
    """Log source enumeration"""
    API = "api"
    DATABASE = "database"
    CELERY = "celery"
    PIPELINE = "pipeline"
    SYSTEM = "system"
    AUTH = "auth"


class LogEntry(BaseModel):
    """Single log entry"""
    timestamp: datetime
    level: LogLevel
    source: LogSource
    message: str
    details: Optional[Dict[str, Any]] = None
    user_id: Optional[str] = None
    ip_address: Optional[str] = None
    request_id: Optional[str] = None


class LogQueryRequest(BaseModel):
    """Query logs request"""
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    levels: Optional[List[LogLevel]] = None
    sources: Optional[List[LogSource]] = None
    user_id: Optional[str] = None
    search: Optional[str] = None  # Search in message
    limit: int = Field(100, ge=1, le=1000)
    offset: int = Field(0, ge=0)


class LogResponse(BaseModel):
    """Log query response"""
    logs: List[LogEntry]
    total: int
    has_more: bool


class LogDownloadRequest(BaseModel):
    """Download logs request"""
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    levels: Optional[List[LogLevel]] = None
    sources: Optional[List[LogSource]] = None
    format: str = "json"  # json, csv, txt


# ============================================================================
# System Health & Maintenance
# ============================================================================

class ServiceStatus(str, Enum):
    """Service status"""
    HEALTHY = "healthy"
    DEGRADED = "degraded"
    DOWN = "down"
    UNKNOWN = "unknown"


class ServiceHealth(BaseModel):
    """Individual service health"""
    name: str
    status: ServiceStatus
    response_time: Optional[float] = None  # in milliseconds
    last_check: datetime
    message: Optional[str] = None
    details: Optional[Dict[str, Any]] = None


class SystemHealth(BaseModel):
    """Overall system health"""
    status: ServiceStatus
    services: List[ServiceHealth]
    timestamp: datetime
    uptime: float  # in seconds
    version: str


class CleanupOptions(BaseModel):
    """System cleanup options"""
    old_logs: bool = False  # Delete logs older than retention period
    temp_files: bool = False  # Delete temporary files
    failed_tasks: bool = False  # Clean up failed task data
    orphaned_files: bool = False  # Remove files not associated with any sample
    old_notifications: bool = False  # Delete old read notifications
    days_to_keep: int = Field(30, ge=1, le=365)
    dry_run: bool = True  # Preview what would be deleted


class CleanupResult(BaseModel):
    """Cleanup operation result"""
    operation: str
    items_deleted: int
    space_freed: int  # bytes
    duration: float  # seconds
    errors: List[str] = []


class CleanupResponse(BaseModel):
    """Complete cleanup response"""
    total_items_deleted: int
    total_space_freed: int
    total_duration: float
    results: List[CleanupResult]
    timestamp: datetime


# ============================================================================
# System Statistics (Admin)
# ============================================================================

class AdminSystemStats(BaseModel):
    """System-wide statistics for admin"""
    total_users: int
    active_users: int
    admin_users: int
    total_projects: int
    total_samples: int
    total_tasks: int
    total_files: int
    total_storage_used: int
    total_storage_allocated: int
    avg_storage_per_user: float
    tasks_today: int
    tasks_this_week: int
    tasks_this_month: int
    success_rate_today: float
    success_rate_this_week: float
    success_rate_this_month: float
    active_sessions: int
    cpu_usage: float
    memory_usage: float
    disk_usage: float
    generated_at: datetime


class UserActivityReport(BaseModel):
    """User activity report"""
    user_id: str
    username: str
    email: str
    last_login: Optional[datetime]
    login_count: int
    projects_created: int
    samples_processed: int
    tasks_run: int
    storage_used: int
    last_activity: Optional[datetime]
    is_active: bool


class SystemActivityReport(BaseModel):
    """System activity report"""
    period: str  # today, week, month, year
    new_users: int
    active_users: int
    projects_created: int
    samples_processed: int
    tasks_executed: int
    tasks_succeeded: int
    tasks_failed: int
    storage_added: int
    peak_concurrent_tasks: int
    avg_task_duration: float
    generated_at: datetime


# ============================================================================
# Audit Log
# ============================================================================

class AuditAction(str, Enum):
    """Audit action types"""
    USER_CREATED = "user_created"
    USER_UPDATED = "user_updated"
    USER_DELETED = "user_deleted"
    USER_ROLE_CHANGED = "user_role_changed"
    USER_ACTIVATED = "user_activated"
    USER_DEACTIVATED = "user_deactivated"
    PASSWORD_RESET = "password_reset"
    CONFIG_UPDATED = "config_updated"
    CONFIG_RESET = "config_reset"
    SYSTEM_CLEANUP = "system_cleanup"
    BACKUP_CREATED = "backup_created"
    BACKUP_RESTORED = "backup_restored"


class AuditLogEntry(BaseModel):
    """Audit log entry"""
    id: str
    timestamp: datetime
    action: AuditAction
    admin_user_id: str
    admin_username: str
    target_user_id: Optional[str] = None
    target_username: Optional[str] = None
    details: Dict[str, Any]
    ip_address: Optional[str] = None
    user_agent: Optional[str] = None


class AuditLogResponse(BaseModel):
    """Audit log query response"""
    logs: List[AuditLogEntry]
    total: int
    page: int
    page_size: int


# ============================================================================
# System Backup
# ============================================================================

class BackupType(str, Enum):
    """Backup type"""
    FULL = "full"
    INCREMENTAL = "incremental"
    DATABASE_ONLY = "database_only"
    FILES_ONLY = "files_only"


class BackupRequest(BaseModel):
    """Create backup request"""
    backup_type: BackupType = BackupType.FULL
    description: Optional[str] = None
    compress: bool = True


class BackupInfo(BaseModel):
    """Backup information"""
    id: str
    backup_type: BackupType
    file_path: str
    size: int  # bytes
    compressed: bool
    description: Optional[str] = None
    created_at: datetime
    created_by: str
    status: str  # completed, in_progress, failed


class BackupListResponse(BaseModel):
    """List of backups"""
    backups: List[BackupInfo]
    total: int


class RestoreRequest(BaseModel):
    """Restore from backup request"""
    backup_id: str
    confirm: bool = True
    restore_database: bool = True
    restore_files: bool = True


# ============================================================================
# Response Messages
# ============================================================================

class AdminOperationResponse(BaseModel):
    """Generic admin operation response"""
    success: bool
    message: str
    details: Optional[Dict[str, Any]] = None
    timestamp: datetime = Field(default_factory=datetime.utcnow)


# ============================================================================
# System Metrics (Enhanced)
# ============================================================================

class SystemMetrics(BaseModel):
    """Detailed system metrics"""
    cpu: Dict[str, Any]  # {"usage": 45.2, "load": [1.2, 1.5, 1.8]}
    memory: Dict[str, Any]  # {"used": 8192, "total": 16384, "usagePercent": 50.0}
    disk: Dict[str, Any]  # {"used": 102400, "total": 512000, "usagePercent": 20.0}
    network: Optional[Dict[str, Any]] = None
    timestamp: datetime = Field(default_factory=datetime.utcnow)


# ============================================================================
# System Alerts
# ============================================================================

class AlertType(str, Enum):
    """Alert type"""
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


class AlertSeverity(str, Enum):
    """Alert severity"""
    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class Alert(BaseModel):
    """System alert"""
    id: str
    type: AlertType
    severity: AlertSeverity
    title: str
    message: str
    timestamp: datetime
    resolved: bool = False
    resolved_at: Optional[datetime] = None
    resolved_by: Optional[str] = None
    source: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class AlertListResponse(BaseModel):
    """Alert list response"""
    alerts: List[Alert]
    total: int
    unresolved_count: int


# ============================================================================
# Resource Management
# ============================================================================

class ResourceUsage(BaseModel):
    """System resource usage"""
    storage: Dict[str, int]  # {"total": 512000, "used": 102400}
    compute: Dict[str, int]  # {"active": 5, "limit": 10}
    memory: Dict[str, int]  # {"total": 16384, "used": 8192}
    timestamp: datetime = Field(default_factory=datetime.utcnow)


# ============================================================================
# Job Management
# ============================================================================

class JobType(str, Enum):
    """Job type"""
    PIPELINE = "pipeline"
    BACKUP = "backup"
    CLEANUP = "cleanup"
    EXPORT = "export"
    IMPORT = "import"
    SYSTEM_TASK = "system_task"


class JobStatus(str, Enum):
    """Job status"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobInfo(BaseModel):
    """Job information"""
    id: str
    type: JobType
    status: JobStatus
    user_id: str
    username: Optional[str] = None
    created_at: datetime
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    progress: Optional[float] = None  # 0.0 to 1.0
    message: Optional[str] = None
    result: Optional[Dict[str, Any]] = None
    error: Optional[str] = None


class JobListResponse(BaseModel):
    """Job list response"""
    jobs: List[JobInfo]
    total: int
    page: int
    page_size: int


class JobOperation(BaseModel):
    """Job operation request"""
    reason: Optional[str] = None


# ============================================================================
# Audit Log Export
# ============================================================================

class AuditLogExportRequest(BaseModel):
    """Audit log export request"""
    start_date: Optional[datetime] = None
    end_date: Optional[datetime] = None
    user_id: Optional[str] = None
    action: Optional[AuditAction] = None
    format: str = "json"  # json, csv


class ExportResult(BaseModel):
    """Export result"""
    download_url: str
    file_name: str
    file_size: int
    expires_at: datetime
