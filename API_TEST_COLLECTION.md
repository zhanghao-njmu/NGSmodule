# NGSmodule API Test Collection

Quick reference for testing backend API endpoints with example requests and expected responses.

## Base Configuration

```bash
BASE_URL="http://localhost:8000/api/v1"
TOKEN="your_jwt_token_here"
```

---

## 1. Authentication

### Register New User

```bash
curl -X POST "${BASE_URL}/auth/register" \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser",
    "email": "test@example.com",
    "password": "Test123!",
    "full_name": "Test User"
  }'
```

**Expected Response** (201):
```json
{
  "id": "uuid",
  "username": "testuser",
  "email": "test@example.com",
  "full_name": "Test User",
  "role": "user",
  "is_active": true,
  "created_at": "2025-11-23T10:00:00Z"
}
```

### Login

```bash
curl -X POST "${BASE_URL}/auth/login" \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser",
    "password": "Test123!"
  }'
```

**Expected Response** (200):
```json
{
  "access_token": "eyJ0eXAiOiJKV1QiLCJhbGc...",
  "token_type": "bearer",
  "user": {
    "id": "uuid",
    "username": "testuser",
    "email": "test@example.com",
    "role": "user"
  }
}
```

**Save token for subsequent requests:**
```bash
export TOKEN="eyJ0eXAiOiJKV1QiLCJhbGc..."
```

### Get Current User

```bash
curl -X GET "${BASE_URL}/auth/me" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Logout

```bash
curl -X POST "${BASE_URL}/auth/logout" \
  -H "Authorization: Bearer ${TOKEN}"
```

---

## 2. Projects (Items)

### List Projects (Paginated)

```bash
curl -X GET "${BASE_URL}/items?skip=0&limit=20" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "items": [
    {
      "id": "uuid",
      "name": "RNA-Seq Project",
      "description": "Differential expression analysis",
      "status": "active",
      "created_at": "2025-11-23T10:00:00Z",
      "sample_count": 10,
      "file_count": 20
    }
  ],
  "total": 1,
  "skip": 0,
  "limit": 20
}
```

### Create Project

```bash
curl -X POST "${BASE_URL}/items" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "New RNA-Seq Project",
    "description": "Testing integration",
    "organism": "Homo sapiens",
    "project_type": "RNA-Seq",
    "status": "active"
  }'
```

**Expected Response** (201):
```json
{
  "id": "uuid",
  "name": "New RNA-Seq Project",
  "description": "Testing integration",
  "organism": "Homo sapiens",
  "project_type": "RNA-Seq",
  "status": "active",
  "created_at": "2025-11-23T10:05:00Z",
  "user_id": "uuid"
}
```

### Get Project by ID

```bash
PROJECT_ID="your_project_id"
curl -X GET "${BASE_URL}/items/${PROJECT_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Update Project

```bash
curl -X PUT "${BASE_URL}/items/${PROJECT_ID}" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Updated Project Name",
    "status": "completed"
  }'
```

### Delete Project

```bash
curl -X DELETE "${BASE_URL}/items/${PROJECT_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "message": "Project deleted successfully",
  "id": "uuid"
}
```

### Get Project Statistics

```bash
curl -X GET "${BASE_URL}/items/${PROJECT_ID}/stats" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "total_samples": 10,
  "total_files": 50,
  "total_tasks": 5,
  "completed_tasks": 3,
  "running_tasks": 1,
  "failed_tasks": 1,
  "storage_used": 1073741824
}
```

---

## 3. Samples

### List Samples

```bash
curl -X GET "${BASE_URL}/samples?project_id=${PROJECT_ID}&skip=0&limit=20" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Create Sample

```bash
curl -X POST "${BASE_URL}/samples" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Sample_001",
    "project_id": "'${PROJECT_ID}'",
    "sample_type": "RNA",
    "organism": "Homo sapiens",
    "tissue": "liver",
    "condition": "control",
    "replicate": 1
  }'
```

**Expected Response** (201):
```json
{
  "id": "uuid",
  "name": "Sample_001",
  "project_id": "uuid",
  "sample_type": "RNA",
  "organism": "Homo sapiens",
  "tissue": "liver",
  "condition": "control",
  "replicate": 1,
  "created_at": "2025-11-23T10:10:00Z"
}
```

### Batch Create Samples

```bash
curl -X POST "${BASE_URL}/samples/batch" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "project_id": "'${PROJECT_ID}'",
    "samples": [
      {
        "name": "Sample_001",
        "sample_type": "RNA",
        "organism": "Homo sapiens",
        "condition": "control"
      },
      {
        "name": "Sample_002",
        "sample_type": "RNA",
        "organism": "Homo sapiens",
        "condition": "treatment"
      }
    ]
  }'
```

**Expected Response** (201):
```json
{
  "created": 2,
  "samples": [...]
}
```

### Update Sample

```bash
SAMPLE_ID="your_sample_id"
curl -X PUT "${BASE_URL}/samples/${SAMPLE_ID}" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "condition": "treatment",
    "replicate": 2
  }'
```

### Delete Sample

```bash
curl -X DELETE "${BASE_URL}/samples/${SAMPLE_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

---

## 4. Files

### List Files

```bash
curl -X GET "${BASE_URL}/files?project_id=${PROJECT_ID}&skip=0&limit=20" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Upload File

```bash
curl -X POST "${BASE_URL}/files/upload" \
  -H "Authorization: Bearer ${TOKEN}" \
  -F "file=@/path/to/sample.fastq.gz" \
  -F "sample_id=${SAMPLE_ID}" \
  -F "file_type=fastq" \
  -F "description=Raw sequencing data"
```

**Expected Response** (201):
```json
{
  "id": "uuid",
  "filename": "sample.fastq.gz",
  "file_type": "fastq",
  "file_size": 1073741824,
  "sample_id": "uuid",
  "storage_path": "s3://bucket/path/sample.fastq.gz",
  "created_at": "2025-11-23T10:15:00Z"
}
```

### Download File

```bash
FILE_ID="your_file_id"
curl -X GET "${BASE_URL}/files/${FILE_ID}/download" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o downloaded_file.fastq.gz
```

### Get File Metadata

```bash
curl -X GET "${BASE_URL}/files/${FILE_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Delete File

```bash
curl -X DELETE "${BASE_URL}/files/${FILE_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

---

## 5. Tasks (Pipeline Execution)

### List Tasks

```bash
curl -X GET "${BASE_URL}/tasks?project_id=${PROJECT_ID}&status=running" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "items": [
    {
      "id": "uuid",
      "task_name": "QC Analysis",
      "task_type": "qc",
      "status": "running",
      "progress": 45,
      "project_id": "uuid",
      "started_at": "2025-11-23T10:20:00Z"
    }
  ],
  "total": 1
}
```

### Create Task

```bash
curl -X POST "${BASE_URL}/tasks" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "task_name": "RNA-Seq Analysis",
    "task_type": "rnaseq",
    "project_id": "'${PROJECT_ID}'",
    "sample_ids": ["'${SAMPLE_ID}'"],
    "parameters": {
      "genome": "hg38",
      "aligner": "STAR",
      "quantifier": "featureCounts"
    }
  }'
```

**Expected Response** (201):
```json
{
  "id": "uuid",
  "task_name": "RNA-Seq Analysis",
  "task_type": "rnaseq",
  "status": "pending",
  "project_id": "uuid",
  "celery_task_id": "celery_uuid",
  "created_at": "2025-11-23T10:25:00Z"
}
```

### Get Task Status

```bash
TASK_ID="your_task_id"
curl -X GET "${BASE_URL}/tasks/${TASK_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "id": "uuid",
  "status": "running",
  "progress": 65,
  "started_at": "2025-11-23T10:20:00Z",
  "logs": [
    {"timestamp": "2025-11-23T10:20:10Z", "message": "Starting analysis..."},
    {"timestamp": "2025-11-23T10:21:00Z", "message": "Quality control complete"}
  ]
}
```

### Cancel Task

```bash
curl -X PUT "${BASE_URL}/tasks/${TASK_ID}/cancel" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Get Task Logs

```bash
curl -X GET "${BASE_URL}/tasks/${TASK_ID}/logs" \
  -H "Authorization: Bearer ${TOKEN}"
```

---

## 6. Results

### List Results

```bash
curl -X GET "${BASE_URL}/results?task_id=${TASK_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "items": [
    {
      "id": "uuid",
      "result_type": "qc_report",
      "task_id": "uuid",
      "file_path": "s3://bucket/results/qc_report.html",
      "metadata": {
        "total_reads": 50000000,
        "quality_score": 35.5
      },
      "created_at": "2025-11-23T10:30:00Z"
    }
  ]
}
```

### Get Result Details

```bash
RESULT_ID="your_result_id"
curl -X GET "${BASE_URL}/results/${RESULT_ID}" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Download Result

```bash
curl -X GET "${BASE_URL}/results/${RESULT_ID}/download" \
  -H "Authorization: Bearer ${TOKEN}" \
  -o result.html
```

---

## 7. Statistics

### Dashboard Statistics

```bash
curl -X GET "${BASE_URL}/stats/dashboard" \
  -H "Authorization: Bearer ${TOKEN}"
```

**Expected Response** (200):
```json
{
  "total_projects": 5,
  "active_projects": 3,
  "total_samples": 50,
  "total_files": 100,
  "running_tasks": 2,
  "completed_tasks": 15,
  "storage_used": 10737418240,
  "storage_quota": 107374182400
}
```

### Project Statistics

```bash
curl -X GET "${BASE_URL}/stats/projects" \
  -H "Authorization: Bearer ${TOKEN}"
```

### Task Statistics

```bash
curl -X GET "${BASE_URL}/stats/tasks" \
  -H "Authorization: Bearer ${TOKEN}"
```

---

## 8. Admin Endpoints (Admin Role Required)

### List All Users

```bash
curl -X GET "${BASE_URL}/admin/users" \
  -H "Authorization: Bearer ${ADMIN_TOKEN}"
```

### Get System Statistics

```bash
curl -X GET "${BASE_URL}/admin/stats" \
  -H "Authorization: Bearer ${ADMIN_TOKEN}"
```

**Expected Response** (200):
```json
{
  "total_users": 25,
  "active_users": 20,
  "total_projects": 100,
  "total_storage_used": 1073741824000,
  "total_storage_quota": 10737418240000,
  "cpu_usage": 45.5,
  "memory_usage": 62.3
}
```

### System Health Check

```bash
curl -X GET "${BASE_URL}/admin/health" \
  -H "Authorization: Bearer ${ADMIN_TOKEN}"
```

**Expected Response** (200):
```json
{
  "overall": "healthy",
  "components": {
    "api": {"status": "operational", "response_time": 25},
    "database": {"status": "operational", "connections": 10},
    "storage": {"status": "operational", "available_space": 500000000000},
    "queue": {"status": "operational", "pending_tasks": 5}
  },
  "metrics": {
    "cpu": 45.5,
    "memory": 62.3,
    "disk": 78.2
  }
}
```

### Update User (Admin)

```bash
USER_ID="target_user_id"
curl -X PUT "${BASE_URL}/admin/users/${USER_ID}" \
  -H "Authorization: Bearer ${ADMIN_TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "role": "admin",
    "storage_quota": 214748364800,
    "is_active": true
  }'
```

### Deactivate User

```bash
curl -X PUT "${BASE_URL}/admin/users/${USER_ID}/deactivate" \
  -H "Authorization: Bearer ${ADMIN_TOKEN}"
```

---

## 9. WebSocket Testing

### Connect to WebSocket

```javascript
// JavaScript example for browser console
const ws = new WebSocket(`ws://localhost:8000/api/v1/ws?token=${TOKEN}`);

ws.onopen = () => {
  console.log('✓ WebSocket connected');
};

ws.onmessage = (event) => {
  const data = JSON.parse(event.data);
  console.log('Received:', data);

  // Expected message format:
  // {
  //   "type": "task_update",
  //   "task_id": "uuid",
  //   "status": "running",
  //   "progress": 75,
  //   "message": "Processing sample 3/4"
  // }
};

ws.onerror = (error) => {
  console.error('WebSocket error:', error);
};

ws.onclose = () => {
  console.log('WebSocket disconnected');
};
```

### Test Task Updates via WebSocket

1. Connect to WebSocket (see above)
2. Create a new task via REST API
3. Observe real-time status updates in WebSocket messages
4. Verify progress updates (0% → 100%)
5. Verify final status message (completed/failed)

---

## 10. Error Response Examples

### 400 Bad Request (Validation Error)

```json
{
  "detail": [
    {
      "loc": ["body", "email"],
      "msg": "value is not a valid email address",
      "type": "value_error.email"
    }
  ]
}
```

### 401 Unauthorized (Missing/Invalid Token)

```json
{
  "detail": "Not authenticated"
}
```

### 403 Forbidden (Insufficient Permissions)

```json
{
  "detail": "Insufficient permissions. Admin role required."
}
```

### 404 Not Found

```json
{
  "detail": "Project not found"
}
```

### 500 Internal Server Error

```json
{
  "detail": "Internal server error. Please contact support."
}
```

---

## Quick Test Script

Save this as `test-api.sh`:

```bash
#!/bin/bash

BASE_URL="http://localhost:8000/api/v1"

# 1. Health check
echo "Testing health endpoint..."
curl -s "${BASE_URL}/health" | jq

# 2. Register user
echo -e "\n\nRegistering new user..."
curl -s -X POST "${BASE_URL}/auth/register" \
  -H "Content-Type: application/json" \
  -d '{
    "username": "testuser'$(date +%s)'",
    "email": "test'$(date +%s)'@example.com",
    "password": "Test123!",
    "full_name": "Test User"
  }' | jq

# 3. Login
echo -e "\n\nLogging in..."
LOGIN_RESPONSE=$(curl -s -X POST "${BASE_URL}/auth/login" \
  -H "Content-Type: application/json" \
  -d '{
    "username": "admin",
    "password": "admin"
  }')

TOKEN=$(echo $LOGIN_RESPONSE | jq -r '.access_token')
echo "Token: ${TOKEN:0:50}..."

# 4. Get current user
echo -e "\n\nGetting current user..."
curl -s "${BASE_URL}/auth/me" \
  -H "Authorization: Bearer ${TOKEN}" | jq

# 5. Create project
echo -e "\n\nCreating project..."
curl -s -X POST "${BASE_URL}/items" \
  -H "Authorization: Bearer ${TOKEN}" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Test Project",
    "description": "Integration test",
    "organism": "Homo sapiens"
  }' | jq

echo -e "\n\nAPI tests complete!"
```

Make executable and run:
```bash
chmod +x test-api.sh
./test-api.sh
```

---

**Note**: Replace placeholder values (UUIDs, tokens) with actual values from your environment.
**Tip**: Use `jq` to format JSON responses for better readability.
