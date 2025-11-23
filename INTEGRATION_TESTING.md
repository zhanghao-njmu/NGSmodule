# NGSmodule Frontend-Backend Integration Testing Guide

## 📋 Overview

This document provides a comprehensive guide for testing the integration between the NGSmodule frontend (React + TypeScript) and backend (FastAPI).

**Status**: ✅ Frontend build successful (0 TypeScript errors)
**Date**: 2025-11-23
**Phase**: 23 Part 5 - Integration Testing

---

## 🏗️ Architecture

### Services
- **Frontend**: React + Vite (Port 3000)
- **Backend API**: FastAPI (Port 8000)
- **Database**: PostgreSQL (Port 5432)
- **Cache**: Redis (Port 6379)
- **Storage**: MinIO (Ports 9000, 9001)
- **Task Queue**: Celery + Flower (Port 5555)

### API Endpoints
- **Base URL**: `http://localhost:8000/api/v1`
- **WebSocket**: `ws://localhost:8000/api/v1/ws`
- **API Docs**: `http://localhost:8000/api/v1/docs`

---

## 🚀 Quick Start

### 1. Environment Setup

```bash
# Frontend environment
cd frontend
cp .env.example .env

# Backend environment
cd ../backend
cp .env.example .env
```

### 2. Start All Services (Docker)

```bash
# From project root
docker-compose up -d

# Check service status
docker-compose ps

# View logs
docker-compose logs -f backend
docker-compose logs -f frontend
```

### 3. Start Development Mode (Local)

**Backend:**
```bash
cd backend

# Install dependencies
pip install -r requirements.txt

# Run migrations
alembic upgrade head

# Create admin user
python create_admin.py

# Start server
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

**Frontend:**
```bash
cd frontend

# Install dependencies
npm install

# Start dev server
npm run dev
```

---

## 🧪 Integration Test Checklist

### Phase 1: Service Health Checks

- [ ] **Backend API** - http://localhost:8000/api/v1/health
- [ ] **Frontend Dev Server** - http://localhost:3000
- [ ] **PostgreSQL** - Connection successful
- [ ] **Redis** - Connection successful
- [ ] **MinIO** - Console accessible at http://localhost:9001
- [ ] **API Documentation** - http://localhost:8000/api/v1/docs

### Phase 2: Authentication Flow

- [ ] **User Registration**
  - [ ] Frontend form validation works
  - [ ] API creates user successfully
  - [ ] Success toast notification appears
  - [ ] Auto-redirect to login

- [ ] **User Login**
  - [ ] Form validation (username + password required)
  - [ ] API returns JWT token
  - [ ] Token stored in localStorage
  - [ ] User redirected to dashboard
  - [ ] authStore updated with user data

- [ ] **Token Refresh**
  - [ ] Token automatically refreshes before expiry
  - [ ] API interceptor adds Authorization header
  - [ ] Failed requests trigger logout on 401

- [ ] **Logout**
  - [ ] Token removed from localStorage
  - [ ] User redirected to login page
  - [ ] authStore reset

### Phase 3: Dashboard & Statistics

- [ ] **Dashboard Load**
  - [ ] Stats API called successfully
  - [ ] Statistics cards display correct data:
    - [ ] Total Projects
    - [ ] Running Tasks
    - [ ] Completed Tasks
    - [ ] Storage Used (%)
  - [ ] Loading skeleton shows during fetch
  - [ ] Error handling if API fails

### Phase 4: Project Management (CRUD)

- [ ] **List Projects**
  - [ ] GET `/items` returns paginated data
  - [ ] Table displays projects correctly
  - [ ] Pagination works (page change triggers API call)
  - [ ] Search filter works
  - [ ] Status filter works
  - [ ] Empty state shows when no projects

- [ ] **Create Project**
  - [ ] Modal opens on "New Project" click
  - [ ] Form validation works (required fields)
  - [ ] POST `/items` creates project
  - [ ] Success toast notification
  - [ ] Table refreshes with new project
  - [ ] Modal closes

- [ ] **Edit Project**
  - [ ] Edit button opens modal with pre-filled data
  - [ ] PUT `/items/{id}` updates project
  - [ ] Table row updates without full refresh
  - [ ] Success notification

- [ ] **Delete Project**
  - [ ] Delete confirmation modal appears
  - [ ] DELETE `/items/{id}` removes project
  - [ ] Table row removed
  - [ ] Success notification

- [ ] **View Project Details**
  - [ ] Click project name navigates to detail page
  - [ ] GET `/items/{id}` fetches full data
  - [ ] All project fields display correctly

### Phase 5: Sample Management

- [ ] **List Samples**
  - [ ] GET `/samples` with project filter
  - [ ] Samples display in grid/list view
  - [ ] Layout toggle works (grid ↔ list)
  - [ ] Group by project dropdown works

- [ ] **Create Sample**
  - [ ] Form validation (sample name, type, organism)
  - [ ] POST `/samples` creates sample
  - [ ] Success notification
  - [ ] Sample appears in list

- [ ] **Batch Import Samples**
  - [ ] CSV/Excel file upload
  - [ ] File parsing and preview
  - [ ] Validation errors highlighted
  - [ ] Bulk create API call
  - [ ] Progress indicator
  - [ ] Summary notification

- [ ] **Edit/Delete Sample**
  - [ ] PUT `/samples/{id}` updates
  - [ ] DELETE `/samples/{id}` removes
  - [ ] Associated files warning

### Phase 6: File Management

- [ ] **List Files**
  - [ ] GET `/files` with filters
  - [ ] Files display with size, type, date
  - [ ] Project/sample filter works
  - [ ] File type icons correct

- [ ] **Upload Files**
  - [ ] File selection (multiple files)
  - [ ] File type validation (.fastq, .bam, etc.)
  - [ ] Upload progress bar
  - [ ] POST `/files/upload` for each file
  - [ ] Success notification per file
  - [ ] List refreshes

- [ ] **Download File**
  - [ ] GET `/files/{id}/download` triggers download
  - [ ] File downloads with correct name
  - [ ] Large file streaming works

- [ ] **Delete File**
  - [ ] Confirmation dialog
  - [ ] DELETE `/files/{id}` removes file
  - [ ] MinIO object deleted
  - [ ] List updates

### Phase 7: Task/Pipeline Execution

- [ ] **List Tasks**
  - [ ] GET `/tasks` returns task list
  - [ ] Status badges correct (pending/running/completed/failed)
  - [ ] Running tasks show progress bar
  - [ ] Refresh button works

- [ ] **Create Task**
  - [ ] Pipeline template selection
  - [ ] Parameter configuration form
  - [ ] Sample/file selection
  - [ ] POST `/tasks` creates task
  - [ ] Task appears in list as "pending"

- [ ] **Real-time Task Updates**
  - [ ] WebSocket connection established
  - [ ] Task status updates automatically
  - [ ] Progress bar updates in real-time
  - [ ] Toast notification on completion/failure
  - [ ] No page refresh needed

- [ ] **View Task Details**
  - [ ] Task logs display
  - [ ] Parameter values shown
  - [ ] Output files listed
  - [ ] Execution timeline

- [ ] **Cancel Task**
  - [ ] PUT `/tasks/{id}/cancel` stops task
  - [ ] Celery task revoked
  - [ ] Status updates to "cancelled"

### Phase 8: Results & Analysis

- [ ] **List Results**
  - [ ] GET `/results` with filters
  - [ ] Results grouped by task
  - [ ] Result type filter (QC/alignment/quantification)
  - [ ] Task ID filter works

- [ ] **View Result**
  - [ ] Result detail page loads
  - [ ] Visualization renders (charts, tables)
  - [ ] Download result button works
  - [ ] Share result link works

- [ ] **Download Results**
  - [ ] GET `/results/{id}/download`
  - [ ] Multiple file formats (JSON, CSV, PDF)
  - [ ] Bulk download for task

### Phase 9: Admin Features

- [ ] **User Management**
  - [ ] GET `/admin/users` lists all users
  - [ ] Admin can view user stats
  - [ ] Edit user (role, quota)
  - [ ] Activate/deactivate user
  - [ ] Delete user (with data cleanup warning)

- [ ] **System Monitoring**
  - [ ] GET `/admin/stats` returns system metrics
  - [ ] Dashboard shows:
    - [ ] CPU/memory usage
    - [ ] Storage usage
    - [ ] Active users
    - [ ] Task queue size
  - [ ] Real-time updates via WebSocket

- [ ] **System Health**
  - [ ] GET `/admin/health` returns component status
  - [ ] Database health
  - [ ] Redis health
  - [ ] MinIO health
  - [ ] Celery worker health

### Phase 10: Error Handling

- [ ] **Network Errors**
  - [ ] API timeout shows error toast
  - [ ] Retry mechanism works
  - [ ] Offline indicator appears

- [ ] **Validation Errors**
  - [ ] Form errors display inline
  - [ ] API validation errors show in toast
  - [ ] Error messages are user-friendly

- [ ] **Authorization Errors**
  - [ ] 401 Unauthorized triggers logout
  - [ ] 403 Forbidden shows access denied message
  - [ ] Protected routes redirect to login

- [ ] **Server Errors**
  - [ ] 500 Internal Server Error shows friendly message
  - [ ] Error details logged to console (dev mode)
  - [ ] Retry option available

### Phase 11: Performance

- [ ] **Loading States**
  - [ ] Skeleton loaders during API calls
  - [ ] Button loading spinners
  - [ ] Table loading overlay

- [ ] **Caching**
  - [ ] API responses cached where appropriate
  - [ ] Cache invalidation on mutations
  - [ ] Stale data refreshed

- [ ] **Pagination**
  - [ ] Large lists paginated (not loading all data)
  - [ ] Virtual scrolling for very long lists
  - [ ] Infinite scroll works smoothly

### Phase 12: UI/UX

- [ ] **Responsive Design**
  - [ ] Mobile layout works (320px+)
  - [ ] Tablet layout works (768px+)
  - [ ] Desktop layout works (1024px+)

- [ ] **Animations**
  - [ ] FadeIn transitions smooth
  - [ ] StaggeredList animations work
  - [ ] Modal open/close animations
  - [ ] No janky animations

- [ ] **Accessibility**
  - [ ] Keyboard navigation works
  - [ ] Focus indicators visible
  - [ ] ARIA labels present
  - [ ] Screen reader friendly

---

## 🐛 Known Issues & Fixes

### Issue Tracking Template

```markdown
## Issue #001: [Brief Description]
**Severity**: Critical | High | Medium | Low
**Component**: Frontend | Backend | Integration
**Status**: Open | In Progress | Fixed

**Description**:
[Detailed description of the issue]

**Steps to Reproduce**:
1.
2.
3.

**Expected Behavior**:
[What should happen]

**Actual Behavior**:
[What actually happens]

**Fix Applied**:
[Description of the fix or N/A if pending]

**Verification**:
- [ ] Issue reproduced
- [ ] Fix implemented
- [ ] Fix tested
- [ ] Issue closed
```

---

## 📊 Testing Report Template

```markdown
# Integration Testing Report
**Date**: YYYY-MM-DD
**Tester**: [Name]
**Environment**: Docker | Local Development
**Frontend**: http://localhost:3000
**Backend**: http://localhost:8000

## Summary
- Total Tests: X
- Passed: X (XX%)
- Failed: X (XX%)
- Blocked: X (XX%)

## Critical Issues
1. [Issue description]

## Non-Critical Issues
1. [Issue description]

## Recommendations
1. [Recommendation]

## Next Steps
1. [Action item]
```

---

## 🔧 Debugging Tips

### Frontend Debugging

```bash
# View API requests in browser console
# Network tab shows all axios requests

# Check store state
console.log(useAuthStore.getState())
console.log(useProjectStore.getState())

# Enable verbose logging
localStorage.setItem('debug', 'app:*')
```

### Backend Debugging

```bash
# View API logs
docker-compose logs -f backend

# Check database
docker-compose exec postgres psql -U ngsmodule -d ngsmodule

# Check Redis
docker-compose exec redis redis-cli

# Check MinIO
# Visit http://localhost:9001 (minioadmin/minioadmin)

# Check Celery tasks
# Visit http://localhost:5555 (Flower)
```

### API Testing with curl

```bash
# Health check
curl http://localhost:8000/api/v1/health

# Login
curl -X POST http://localhost:8000/api/v1/auth/login \
  -H "Content-Type: application/json" \
  -d '{"username":"admin","password":"admin"}'

# Get projects (with token)
curl http://localhost:8000/api/v1/items \
  -H "Authorization: Bearer YOUR_TOKEN"
```

---

## ✅ Sign-off

Once all tests pass:

- [ ] All critical functionality works
- [ ] No console errors
- [ ] No network errors
- [ ] Performance acceptable
- [ ] UI/UX polished
- [ ] Documentation updated
- [ ] Ready for production deployment

**Approved by**: _________________
**Date**: _________________

---

## 📚 Additional Resources

- [FastAPI Documentation](https://fastapi.tiangolo.com/)
- [React Testing Library](https://testing-library.com/react)
- [Ant Design Components](https://ant.design/components/)
- [Zustand State Management](https://github.com/pmndrs/zustand)
- [Axios HTTP Client](https://axios-http.com/)

---

**Generated**: Phase 23 Part 5A - Frontend-Backend Integration Testing
**Next**: Execute tests and document results
