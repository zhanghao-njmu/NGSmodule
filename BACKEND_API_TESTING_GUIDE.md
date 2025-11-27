# 🧪 Backend API Testing & Deployment Guide

## 📋 Table of Contents

1. [Environment Setup](#environment-setup)
2. [Database Initialization](#database-initialization)
3. [Running the Backend](#running-the-backend)
4. [API Testing](#api-testing)
5. [Stats API Endpoints](#stats-api-endpoints)
6. [Notifications API Endpoints](#notifications-api-endpoints)
7. [Admin API Endpoints](#admin-api-endpoints)
8. [Production Deployment](#production-deployment)
9. [Troubleshooting](#troubleshooting)

---

## 🔧 Environment Setup

### Prerequisites

- Python 3.11+
- PostgreSQL 14+
- Redis 7+
- MinIO (for object storage)

### 1. Install Dependencies

```bash
cd backend

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### 2. Configure Environment Variables

```bash
# Copy environment template
cp .env.example .env

# Edit .env with your configuration
nano .env
```

Key variables to configure:
```env
# Database
POSTGRES_USER=ngsmodule
POSTGRES_PASSWORD=your_secure_password
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_DB=ngsmodule

# Redis
REDIS_HOST=localhost
REDIS_PORT=6379
REDIS_PASSWORD=your_redis_password

# Security
SECRET_KEY=your_secret_key_here
JWT_SECRET_KEY=your_jwt_secret_here

# MinIO
MINIO_ROOT_USER=minioadmin
MINIO_ROOT_PASSWORD=your_minio_password
```

---

## 🗄️ Database Initialization

### Option 1: Fresh Installation (New Database)

If you're setting up a fresh database:

```bash
cd backend

# Initialize database with all tables
python init_db.py
```

This will create all tables including:
- users
- projects
- samples
- files
- pipeline_tasks
- results
- pipeline_templates
- **notifications** (new)
- **notification_settings** (new)

### Option 2: Existing Database (Migration)

If you have an existing database and want to add the new notification tables:

```bash
cd backend

# Run Alembic migration
alembic upgrade head
```

This will:
- Add `notifications` table
- Add `notification_settings` table
- Add `last_login` column to `users` table

### Verify Database Setup

```bash
# Connect to PostgreSQL
psql -U ngsmodule -d ngsmodule -h localhost

# List all tables
\dt

# Should see:
# - users
# - projects
# - samples
# - files
# - pipeline_tasks
# - results
# - pipeline_templates
# - notifications
# - notification_settings
# - alembic_version

# Exit psql
\q
```

---

## 🚀 Running the Backend

### Development Mode

```bash
cd backend

# Activate virtual environment
source venv/bin/activate

# Run with auto-reload
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

### Production Mode

```bash
cd backend

# Run with multiple workers
uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4
```

### Verify Backend is Running

```bash
# Check health endpoint
curl http://localhost:8000/health

# Expected response:
# {"status": "healthy", "timestamp": "2024-01-01T12:00:00"}

# Check API docs
open http://localhost:8000/docs
```

---

## 🧪 API Testing

### Authentication

Most endpoints require authentication. First, obtain an access token:

#### 1. Register a User

```bash
curl -X POST http://localhost:8000/api/v1/auth/register \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser",
    "email": "test@example.com",
    "password": "SecurePassword123!",
    "full_name": "Test User"
  }'
```

#### 2. Login

```bash
curl -X POST http://localhost:8000/api/v1/auth/login \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser",
    "password": "SecurePassword123!"
  }'
```

Response:
```json
{
  "access_token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "refresh_token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "token_type": "bearer"
}
```

Save the `access_token` for subsequent requests:
```bash
export TOKEN="eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9..."
```

---

## 📊 Stats API Endpoints

### 1. Get Summary Statistics

Get comprehensive statistics across all entities.

```bash
curl -X GET http://localhost:8000/api/v1/stats/summary \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "projects": {
    "total": 10,
    "active": 7,
    "completed": 2,
    "failed": 1
  },
  "samples": {
    "total": 45,
    "processed": 30,
    "processing": 10,
    "failed": 5
  },
  "tasks": {
    "total": 120,
    "pending": 20,
    "running": 15,
    "completed": 80,
    "failed": 5,
    "success_rate": 94.12
  },
  "files": {
    "total": 200,
    "size": 524288000000,
    "by_type": {
      "fastq": 80,
      "bam": 40,
      "vcf": 30,
      "other": 50
    }
  },
  "storage": {
    "total": 107374182400,
    "used": 52428800000,
    "available": 54945382400,
    "percent_used": 48.8
  },
  "generated_at": "2024-01-01T12:00:00Z"
}
```

### 2. Get Project Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/projects \
  -H "Authorization: Bearer $TOKEN"
```

### 3. Get Quick Stats (Dashboard)

Optimized endpoint for dashboard display:

```bash
curl -X GET http://localhost:8000/api/v1/stats/quick \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "total_projects": 10,
  "active_projects": 7,
  "total_samples": 45,
  "total_tasks": 120,
  "running_tasks": 15,
  "completed_tasks": 80,
  "storage_used": 52428800000,
  "storage_quota": 107374182400,
  "storage_percent": 48.8,
  "unread_notifications": 5
}
```

### 4. Get Trend Data

Get time-series data for metrics:

```bash
# Get task trends for last 30 days
curl -X GET "http://localhost:8000/api/v1/stats/trends/tasks?period=daily&days=30" \
  -H "Authorization: Bearer $TOKEN"

# Get storage trends for last 7 days
curl -X GET "http://localhost:8000/api/v1/stats/trends/storage?period=daily&days=7" \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "metric": "tasks",
  "period": "daily",
  "data": [
    {
      "timestamp": "2024-01-01T00:00:00Z",
      "value": 15,
      "label": "2024-01-01"
    },
    {
      "timestamp": "2024-01-02T00:00:00Z",
      "value": 20,
      "label": "2024-01-02"
    }
  ],
  "total": 2
}
```

### 5. Get Storage Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/storage \
  -H "Authorization: Bearer $TOKEN"
```

### 6. Get Sample Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/samples \
  -H "Authorization: Bearer $TOKEN"
```

### 7. Get Task Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/tasks \
  -H "Authorization: Bearer $TOKEN"
```

### 8. Get File Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/files \
  -H "Authorization: Bearer $TOKEN"
```

### 9. Get Pipeline Statistics

```bash
curl -X GET http://localhost:8000/api/v1/stats/pipelines \
  -H "Authorization: Bearer $TOKEN"
```

### 10. Get System Statistics (Admin Only)

```bash
curl -X GET http://localhost:8000/api/v1/stats/system \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "cpu_usage": 45.2,
  "memory_usage": 68.5,
  "disk_usage": 72.3,
  "active_connections": 25,
  "uptime_seconds": 86400
}
```

### 11. Get User Activity Stats (Admin Only)

```bash
curl -X GET http://localhost:8000/api/v1/stats/users \
  -H "Authorization: Bearer $TOKEN"
```

---

## 🔔 Notifications API Endpoints

### 1. Get Notifications (Paginated)

```bash
# Get all notifications
curl -X GET "http://localhost:8000/api/v1/notifications?skip=0&limit=20" \
  -H "Authorization: Bearer $TOKEN"

# Get only unread notifications
curl -X GET "http://localhost:8000/api/v1/notifications?unread_only=true" \
  -H "Authorization: Bearer $TOKEN"

# Filter by type
curl -X GET "http://localhost:8000/api/v1/notifications?type=task_completed" \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "items": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440000",
      "user_id": "550e8400-e29b-41d4-a716-446655440001",
      "type": "task_completed",
      "title": "Task Completed",
      "message": "Your RNA-seq analysis has completed successfully",
      "data": {
        "task_id": "550e8400-e29b-41d4-a716-446655440002",
        "task_name": "RNA-seq Analysis",
        "status": "completed"
      },
      "read": false,
      "action_url": "/tasks/550e8400-e29b-41d4-a716-446655440002",
      "priority": "normal",
      "created_at": "2024-01-01T12:00:00Z",
      "read_at": null
    }
  ],
  "total": 15,
  "page": 1,
  "page_size": 20,
  "unread_count": 5
}
```

### 2. Get Unread Count

```bash
curl -X GET http://localhost:8000/api/v1/notifications/unread/count \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "count": 5,
  "by_type": {
    "task_completed": 2,
    "task_failed": 1,
    "system_alert": 2
  }
}
```

### 3. Get Single Notification

```bash
curl -X GET http://localhost:8000/api/v1/notifications/550e8400-e29b-41d4-a716-446655440000 \
  -H "Authorization: Bearer $TOKEN"
```

### 4. Mark as Read

```bash
curl -X PUT http://localhost:8000/api/v1/notifications/550e8400-e29b-41d4-a716-446655440000/read \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "id": "550e8400-e29b-41d4-a716-446655440000",
  "read": true,
  "read_at": "2024-01-01T12:30:00Z"
}
```

### 5. Mark All as Read

```bash
curl -X PUT http://localhost:8000/api/v1/notifications/read-all \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "marked_count": 5,
  "message": "All notifications marked as read"
}
```

### 6. Delete Notification

```bash
curl -X DELETE http://localhost:8000/api/v1/notifications/550e8400-e29b-41d4-a716-446655440000 \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "message": "Notification deleted successfully"
}
```

### 7. Get Notification Settings

```bash
curl -X GET http://localhost:8000/api/v1/notifications/settings/current \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "id": "550e8400-e29b-41d4-a716-446655440003",
  "user_id": "550e8400-e29b-41d4-a716-446655440001",
  "email_enabled": true,
  "email_task_completed": true,
  "email_task_failed": true,
  "email_system_alerts": true,
  "app_enabled": true,
  "app_task_updates": true,
  "app_project_updates": true,
  "app_system_alerts": true,
  "push_enabled": false,
  "updated_at": "2024-01-01T10:00:00Z"
}
```

### 8. Update Notification Settings

```bash
curl -X PUT http://localhost:8000/api/v1/notifications/settings/current \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "email_enabled": true,
    "email_task_completed": true,
    "email_task_failed": true,
    "app_enabled": true,
    "app_task_updates": false,
    "push_enabled": false
  }'
```

### 9. Create Notification (Programmatic)

This is typically called from backend services, not directly by users:

```python
# Python example
from app.services.notification_service import NotificationService
from app.schemas.notification import NotificationCreate

notification_service = NotificationService(db)

# Create task notification
notification_service.create_task_notification(
    user_id=user.id,
    task_id=task.id,
    task_name=task.name,
    status="completed"
)

# Create system notification
notification_service.create_system_notification(
    user_id=user.id,
    title="System Maintenance",
    message="Scheduled maintenance will occur tonight at 10 PM",
    priority="high"
)
```

---

## 👨‍💼 Admin API Endpoints

### Admin Authentication

Admin endpoints require admin role. First, create an admin user or promote existing user to admin.

#### Promote User to Admin (via database)

```bash
psql -U ngsmodule -d ngsmodule -c "UPDATE users SET role='admin' WHERE username='testuser';"
```

### 1. Get All Users (Paginated)

```bash
curl -X GET "http://localhost:8000/api/v1/admin/users?skip=0&limit=50" \
  -H "Authorization: Bearer $TOKEN"
```

**Query Parameters:**
- `skip`: Pagination offset (default: 0)
- `limit`: Page size (default: 50, max: 100)
- `role`: Filter by role (user, admin)
- `is_active`: Filter by status (true, false)
- `search`: Search in username, email, name, org
- `sort_by`: Sort field (created_at, username, email, role, storage_used)
- `sort_order`: Sort order (asc, desc)

**Response:**
```json
{
  "users": [
    {
      "id": "550e8400-e29b-41d4-a716-446655440001",
      "username": "testuser",
      "email": "test@example.com",
      "full_name": "Test User",
      "role": "user",
      "organization": "Test Org",
      "is_active": true,
      "storage_used": 52428800,
      "storage_quota": 107374182400,
      "storage_percent": 0.05,
      "last_login": "2024-01-01T12:00:00Z",
      "created_at": "2024-01-01T10:00:00Z",
      "updated_at": "2024-01-01T12:00:00Z"
    }
  ],
  "total": 10,
  "page": 1,
  "page_size": 50,
  "total_pages": 1
}
```

### 2. Get User Details

```bash
curl -X GET http://localhost:8000/api/v1/admin/users/USER_ID \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "id": "550e8400-e29b-41d4-a716-446655440001",
  "username": "testuser",
  "email": "test@example.com",
  "total_projects": 5,
  "total_samples": 20,
  "total_tasks": 50,
  "total_files": 100,
  "last_activity": "2024-01-01T15:00:00Z",
  "login_count": 25
}
```

### 3. Update User Information

```bash
curl -X PUT http://localhost:8000/api/v1/admin/users/USER_ID \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "full_name": "Updated Name",
    "storage_quota": 214748364800,
    "organization": "New Org"
  }'
```

### 4. Change User Role

```bash
# Promote to admin
curl -X PUT http://localhost:8000/api/v1/admin/users/USER_ID/role \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "role": "admin"
  }'

# Demote to user
curl -X PUT http://localhost:8000/api/v1/admin/users/USER_ID/role \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "role": "user"
  }'
```

**Note:** Admins cannot change their own role.

### 5. Activate/Deactivate User

```bash
# Deactivate user
curl -X PUT http://localhost:8000/api/v1/admin/users/USER_ID/activate \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "is_active": false,
    "reason": "Account suspended for policy violation"
  }'

# Activate user
curl -X PUT http://localhost:8000/api/v1/admin/users/USER_ID/activate \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "is_active": true,
    "reason": "Reinstatement after review"
  }'
```

**Note:** Admins cannot deactivate themselves.

### 6. Reset User Password

```bash
curl -X POST http://localhost:8000/api/v1/admin/users/USER_ID/reset-password \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "new_password": "NewSecure123!",
    "notify_user": true
  }'
```

**Response:**
```json
{
  "success": true,
  "message": "Password reset successfully",
  "details": {
    "notify_user": true
  },
  "timestamp": "2024-01-01T12:00:00Z"
}
```

### 7. Delete User

```bash
curl -X DELETE http://localhost:8000/api/v1/admin/users/USER_ID \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "confirm": true,
    "transfer_data_to": "OTHER_USER_ID",
    "reason": "User requested account deletion"
  }'
```

**Warning:** This action is irreversible! Use `transfer_data_to` to transfer user's data to another user before deletion.

### 8. Get System Configuration

```bash
curl -X GET http://localhost:8000/api/v1/admin/config \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "general": {
    "app_name": "NGSmodule",
    "app_version": "1.0.0",
    "maintenance_mode": false,
    "allow_registration": true,
    "default_user_quota": 107374182400,
    "max_upload_size": 524288000,
    "session_timeout": 3600
  },
  "security": {
    "jwt_expiry": 1800,
    "refresh_token_expiry": 604800,
    "password_min_length": 8,
    "require_email_verification": false,
    "enable_2fa": false,
    "max_login_attempts": 5,
    "lockout_duration": 900
  },
  "storage": {
    "storage_backend": "minio",
    "storage_path": "/data/storage",
    "enable_compression": true,
    "retention_policy_days": 365
  },
  "last_updated": "2024-01-01T12:00:00Z"
}
```

### 9. Update System Configuration

```bash
curl -X PUT http://localhost:8000/api/v1/admin/config \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "category": "general",
    "updates": {
      "default_user_quota": 214748364800,
      "max_upload_size": 1048576000
    },
    "reason": "Increased user quotas and upload limits"
  }'
```

### 10. Reset System Configuration

```bash
# Reset all configuration
curl -X POST http://localhost:8000/api/v1/admin/config/reset \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "confirm": true
  }'

# Reset specific categories
curl -X POST http://localhost:8000/api/v1/admin/config/reset \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "categories": ["security", "performance"],
    "confirm": true
  }'
```

**Warning:** This will reset configuration to factory defaults!

### 11. Query System Logs

```bash
# Get recent logs
curl -X GET "http://localhost:8000/api/v1/admin/logs?limit=100" \
  -H "Authorization: Bearer $TOKEN"

# Filter by level
curl -X GET "http://localhost:8000/api/v1/admin/logs?levels=error&levels=critical" \
  -H "Authorization: Bearer $TOKEN"

# Filter by source
curl -X GET "http://localhost:8000/api/v1/admin/logs?sources=api&sources=database" \
  -H "Authorization: Bearer $TOKEN"

# Search logs
curl -X GET "http://localhost:8000/api/v1/admin/logs?search=authentication&limit=50" \
  -H "Authorization: Bearer $TOKEN"

# Date range
curl -X GET "http://localhost:8000/api/v1/admin/logs?start_date=2024-01-01T00:00:00Z&end_date=2024-01-02T00:00:00Z" \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "logs": [
    {
      "timestamp": "2024-01-01T12:00:00Z",
      "level": "error",
      "source": "api",
      "message": "Failed to connect to external service",
      "details": {
        "request_id": "req_123",
        "duration_ms": 150
      },
      "user_id": null,
      "ip_address": "192.168.1.1"
    }
  ],
  "total": 1,
  "has_more": false
}
```

### 12. Download System Logs

```bash
# Download as JSON
curl -X GET "http://localhost:8000/api/v1/admin/logs/download?format=json" \
  -H "Authorization: Bearer $TOKEN" \
  -o logs.json

# Download as CSV
curl -X GET "http://localhost:8000/api/v1/admin/logs/download?format=csv" \
  -H "Authorization: Bearer $TOKEN" \
  -o logs.csv

# Download as text
curl -X GET "http://localhost:8000/api/v1/admin/logs/download?format=txt" \
  -H "Authorization: Bearer $TOKEN" \
  -o logs.txt

# Download with filters
curl -X GET "http://localhost:8000/api/v1/admin/logs/download?format=json&levels=error&start_date=2024-01-01T00:00:00Z" \
  -H "Authorization: Bearer $TOKEN" \
  -o error_logs.json
```

### 13. Get System Health

```bash
curl -X GET http://localhost:8000/api/v1/admin/system/health \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "status": "healthy",
  "services": [
    {
      "name": "PostgreSQL",
      "status": "healthy",
      "response_time": 5.2,
      "last_check": "2024-01-01T12:00:00Z",
      "message": "Database is responsive"
    },
    {
      "name": "Redis",
      "status": "healthy",
      "response_time": 1.5,
      "last_check": "2024-01-01T12:00:00Z",
      "message": "Cache is operational"
    },
    {
      "name": "MinIO",
      "status": "healthy",
      "response_time": 8.3,
      "last_check": "2024-01-01T12:00:00Z",
      "message": "Object storage is available"
    },
    {
      "name": "Celery",
      "status": "healthy",
      "response_time": 3.1,
      "last_check": "2024-01-01T12:00:00Z",
      "message": "Task queue is processing"
    }
  ],
  "timestamp": "2024-01-01T12:00:00Z",
  "uptime": 86400.0,
  "version": "1.0.0"
}
```

### 14. System Cleanup

```bash
# Dry run (preview what would be deleted)
curl -X POST http://localhost:8000/api/v1/admin/system/cleanup \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "old_logs": true,
    "temp_files": true,
    "failed_tasks": true,
    "old_notifications": true,
    "days_to_keep": 30,
    "dry_run": true
  }'

# Actual cleanup
curl -X POST http://localhost:8000/api/v1/admin/system/cleanup \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "old_logs": true,
    "temp_files": true,
    "failed_tasks": true,
    "orphaned_files": false,
    "old_notifications": true,
    "days_to_keep": 30,
    "dry_run": false
  }'
```

**Response:**
```json
{
  "total_items_deleted": 150,
  "total_space_freed": 524288000,
  "total_duration": 2.5,
  "results": [
    {
      "operation": "old_logs",
      "items_deleted": 50,
      "space_freed": 104857600,
      "duration": 0.5,
      "errors": []
    },
    {
      "operation": "temp_files",
      "items_deleted": 30,
      "space_freed": 209715200,
      "duration": 0.3,
      "errors": []
    },
    {
      "operation": "failed_tasks",
      "items_deleted": 20,
      "space_freed": 0,
      "duration": 0.8,
      "errors": []
    },
    {
      "operation": "old_notifications",
      "items_deleted": 50,
      "space_freed": 0,
      "duration": 0.4,
      "errors": []
    }
  ],
  "timestamp": "2024-01-01T12:00:00Z"
}
```

**Recommendation:** Always run with `dry_run: true` first to preview changes!

### 15. Get System Statistics (Admin)

```bash
curl -X GET http://localhost:8000/api/v1/admin/system/stats \
  -H "Authorization: Bearer $TOKEN"
```

**Response:**
```json
{
  "total_users": 100,
  "active_users": 85,
  "admin_users": 5,
  "total_projects": 250,
  "total_samples": 1500,
  "total_tasks": 5000,
  "total_files": 10000,
  "total_storage_used": 5497558138880,
  "total_storage_allocated": 10737418240000,
  "avg_storage_per_user": 54975581388,
  "tasks_today": 50,
  "tasks_this_week": 350,
  "tasks_this_month": 1200,
  "success_rate_today": 92.5,
  "success_rate_this_week": 94.2,
  "success_rate_this_month": 93.8,
  "active_sessions": 15,
  "cpu_usage": 45.2,
  "memory_usage": 68.5,
  "disk_usage": 72.3,
  "generated_at": "2024-01-01T12:00:00Z"
}
```

---

## 🐳 Production Deployment

### Using Docker (Recommended)

#### 1. Setup Production Environment

```bash
# Run the setup script
./deploy-production.sh setup

# This will:
# - Create .env file from .env.production template
# - Generate secure random secrets
# - Create required directories (backups, ssl, logs)
# - Set proper permissions
```

#### 2. Update Production Configuration

Edit `.env` file and update:
- Domain names (replace `yourdomain.com`)
- Email settings (if using notifications)
- SSL certificate paths
- Other production-specific settings

#### 3. Setup SSL Certificates

**Option A: Let's Encrypt (Recommended)**
```bash
certbot certonly --standalone -d yourdomain.com
cp /etc/letsencrypt/live/yourdomain.com/fullchain.pem ssl/cert.pem
cp /etc/letsencrypt/live/yourdomain.com/privkey.pem ssl/key.pem
```

**Option B: Self-Signed Certificate (Testing)**
```bash
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
  -keyout ssl/key.pem -out ssl/cert.pem
```

#### 4. Build Docker Images

```bash
./deploy-production.sh build
```

#### 5. Start Services

```bash
./deploy-production.sh start
```

This starts:
- PostgreSQL database
- Redis cache
- MinIO object storage
- Backend API (with multiple workers)
- Frontend (Nginx)
- Celery workers
- Flower (Celery monitoring)

#### 6. Initialize Database

**For new deployment:**
```bash
# Access backend container
docker-compose -f docker-compose.prod.yml exec backend bash

# Initialize database
python init_db.py

# Exit container
exit
```

**For existing deployment (migration):**
```bash
# Run migration
docker-compose -f docker-compose.prod.yml exec backend alembic upgrade head
```

#### 7. Verify Deployment

```bash
# Check service status
./deploy-production.sh status

# View logs
./deploy-production.sh logs

# Check health endpoints
curl https://yourdomain.com/api/health
curl https://yourdomain.com/api/v1/stats/quick
```

### Manual Deployment

If not using Docker:

#### 1. Setup Services

```bash
# Install PostgreSQL
sudo apt-get install postgresql-14

# Install Redis
sudo apt-get install redis-server

# Install MinIO
wget https://dl.min.io/server/minio/release/linux-amd64/minio
chmod +x minio
sudo mv minio /usr/local/bin/
```

#### 2. Configure Services

```bash
# PostgreSQL
sudo -u postgres createuser ngsmodule
sudo -u postgres createdb ngsmodule
sudo -u postgres psql -c "ALTER USER ngsmodule WITH PASSWORD 'your_password';"

# Redis
sudo nano /etc/redis/redis.conf
# Set: requirepass your_redis_password

# MinIO
minio server /data/minio --console-address ":9001" &
```

#### 3. Setup Backend

```bash
cd backend

# Install dependencies
pip install -r requirements.txt

# Initialize database
python init_db.py

# Run with systemd
sudo nano /etc/systemd/system/ngsmodule-backend.service
```

**systemd service file:**
```ini
[Unit]
Description=NGSmodule Backend API
After=network.target postgresql.service redis.service

[Service]
Type=notify
User=www-data
Group=www-data
WorkingDirectory=/opt/ngsmodule/backend
Environment="PATH=/opt/ngsmodule/backend/venv/bin"
ExecStart=/opt/ngsmodule/backend/venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4
Restart=always

[Install]
WantedBy=multi-user.target
```

```bash
# Enable and start service
sudo systemctl enable ngsmodule-backend
sudo systemctl start ngsmodule-backend
```

#### 4. Setup Nginx

```bash
sudo nano /etc/nginx/sites-available/ngsmodule
```

**Nginx configuration:**
```nginx
server {
    listen 80;
    server_name yourdomain.com;
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl http2;
    server_name yourdomain.com;

    ssl_certificate /path/to/ssl/cert.pem;
    ssl_certificate_key /path/to/ssl/key.pem;

    # Frontend
    location / {
        root /opt/ngsmodule/frontend/dist;
        try_files $uri $uri/ /index.html;
    }

    # Backend API
    location /api {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # WebSocket
    location /ws {
        proxy_pass http://localhost:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

```bash
# Enable site
sudo ln -s /etc/nginx/sites-available/ngsmodule /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

---

## 🛠️ Troubleshooting

### Database Connection Issues

```bash
# Check if PostgreSQL is running
sudo systemctl status postgresql

# Test connection
psql -U ngsmodule -d ngsmodule -h localhost

# Check logs
sudo tail -f /var/log/postgresql/postgresql-14-main.log
```

### Backend Not Starting

```bash
# Check logs
docker-compose -f docker-compose.prod.yml logs backend

# Or if running manually:
tail -f /opt/ngsmodule/backend/logs/app.log

# Common issues:
# 1. Database connection refused - Check POSTGRES_HOST and credentials
# 2. Redis connection refused - Check REDIS_HOST and password
# 3. Port already in use - Check if port 8000 is available
```

### API Returns 500 Error

```bash
# Check backend logs for stack trace
docker-compose -f docker-compose.prod.yml logs backend | grep ERROR

# Common issues:
# 1. Missing environment variables
# 2. Database table not created
# 3. Invalid JWT token
```

### Migration Fails

```bash
# Check current revision
alembic current

# Check migration history
alembic history

# Reset to base (WARNING: destroys data)
alembic downgrade base

# Apply all migrations
alembic upgrade head

# If migration fails, check:
# 1. Database connection
# 2. Table already exists (use alembic stamp head)
# 3. Missing dependencies
```

### Stats API Returns Empty Data

```bash
# Ensure database has data
psql -U ngsmodule -d ngsmodule -c "SELECT COUNT(*) FROM projects;"
psql -U ngsmodule -d ngsmodule -c "SELECT COUNT(*) FROM samples;"
psql -U ngsmodule -d ngsmodule -c "SELECT COUNT(*) FROM pipeline_tasks;"

# If empty, create sample data (development only)
python scripts/seed_database.py
```

### Notifications Not Appearing

```bash
# Check if notifications table exists
psql -U ngsmodule -d ngsmodule -c "\d notifications"

# Check notification settings
curl -X GET http://localhost:8000/api/v1/notifications/settings/current \
  -H "Authorization: Bearer $TOKEN"

# Verify user has notification settings
psql -U ngsmodule -d ngsmodule -c "SELECT * FROM notification_settings;"
```

### Performance Issues

```bash
# Check database query performance
psql -U ngsmodule -d ngsmodule
# Enable query logging
SET log_statement = 'all';

# Check connection pool
# Add to main.py:
# from app.core.database import engine
# print(f"Pool size: {engine.pool.size()}")
# print(f"Pool checked out: {engine.pool.checkedout()}")

# Increase worker count
# Edit docker-compose.prod.yml:
# environment:
#   UVICORN_WORKERS: 8  # Increase based on CPU cores
```

---

## 📚 Additional Resources

### API Documentation

- **Interactive API Docs**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

### Monitoring

- **Flower (Celery)**: http://localhost:5555
- **MinIO Console**: http://localhost:9001

### Database Management

```bash
# Backup database
./deploy-production.sh backup

# Restore database
./deploy-production.sh restore backups/postgres/backup_20240101_120000.sql.gz

# View backups
ls -lh backups/postgres/
```

### Logs

```bash
# View all logs
./deploy-production.sh logs

# View specific service logs
docker-compose -f docker-compose.prod.yml logs backend
docker-compose -f docker-compose.prod.yml logs postgres
docker-compose -f docker-compose.prod.yml logs redis

# Follow logs in real-time
docker-compose -f docker-compose.prod.yml logs -f backend
```

---

## ✅ Testing Checklist

Before deploying to production:

- [ ] Database initialization completed successfully
- [ ] All migrations applied (`alembic upgrade head`)
- [ ] Backend starts without errors
- [ ] Health endpoint responds: `/health`
- [ ] Authentication works (register, login)
- [ ] Stats API returns data: `/api/v1/stats/summary`
- [ ] Notifications API works: `/api/v1/notifications`
- [ ] Frontend can connect to backend
- [ ] WebSocket connection established
- [ ] File upload works
- [ ] Pipeline execution works
- [ ] SSL certificates configured
- [ ] Environment variables set correctly
- [ ] Firewall rules configured
- [ ] Backup strategy in place
- [ ] Monitoring configured

---

## 🎯 Quick Test Script

```bash
#!/bin/bash
# Test all API endpoints quickly

export API_URL="http://localhost:8000"
export TOKEN="" # Add your token here

# Test authentication
echo "Testing authentication..."
curl -X POST $API_URL/api/v1/auth/login \
  -H "Content-Type: application/json" \
  -d '{"username":"testuser","password":"password"}' | jq

# Test stats
echo -e "\n\nTesting stats API..."
curl -X GET $API_URL/api/v1/stats/summary \
  -H "Authorization: Bearer $TOKEN" | jq

curl -X GET $API_URL/api/v1/stats/quick \
  -H "Authorization: Bearer $TOKEN" | jq

# Test notifications
echo -e "\n\nTesting notifications API..."
curl -X GET $API_URL/api/v1/notifications \
  -H "Authorization: Bearer $TOKEN" | jq

curl -X GET $API_URL/api/v1/notifications/unread/count \
  -H "Authorization: Bearer $TOKEN" | jq

curl -X GET $API_URL/api/v1/notifications/settings/current \
  -H "Authorization: Bearer $TOKEN" | jq

echo -e "\n\n✅ All tests completed!"
```

Save as `test-api.sh`, make executable, and run:
```bash
chmod +x test-api.sh
./test-api.sh
```

---

## 📞 Support

For issues or questions:
1. Check logs first: `./deploy-production.sh logs`
2. Review this guide's troubleshooting section
3. Check API documentation: http://localhost:8000/docs
4. Contact development team

---

**Last Updated**: 2024-01-01
**Version**: 1.0.0
