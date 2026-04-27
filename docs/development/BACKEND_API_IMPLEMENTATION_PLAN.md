# Backend API 实现计划

**项目**: NGSmodule 后端API实现
**目标**: 完成前后端API对接，实现所有前端需求的端点
**日期**: 2024-01-01

---

## 📊 现状分析

### 后端已实现的API模块 ✅

**基础模块** (9个):
1. ✅ `/api/v1/auth` - 认证 (登录、注册、token)
2. ✅ `/api/v1/users` - 用户管理
3. ✅ `/api/v1/projects` - 项目管理
4. ✅ `/api/v1/samples` - 样本管理
5. ✅ `/api/v1/files` - 文件管理
6. ✅ `/api/v1/tasks` - 任务管理
7. ✅ `/api/v1/pipelines` - 流程管理
8. ✅ `/api/v1/results` - 结果管理
9. ✅ `/api/v1/ws` - WebSocket实时通信

### 前端需要但后端缺失的API模块 ⚠️

**缺失模块** (4个):
1. ❌ `/api/v1/ai` - AI智能功能 (23个方法)
2. ❌ `/api/v1/analytics` - 数据分析统计
3. ❌ `/api/v1/notifications` - 通知管理
4. ❌ `/api/v1/admin` - 管理后台功能
5. ❌ `/api/v1/stats` - 统计数据

---

## 🎯 前端API需求分析

### 1. AI Intelligence API (`ai.service.ts`)

**方法数**: 23个方法
**复杂度**: 高
**优先级**: P1 (核心功能)

#### A. 参数推荐 (3个端点)

```typescript
// POST /api/v1/ai/recommendations/parameters
getParameterRecommendations(request: RecommendationRequest)
  → PipelineRecommendation

// POST /api/v1/ai/recommendations/parameters/{pipelineType}/{parameterName}
getParameterOptions(pipelineType, parameterName, context)
  → ParameterRecommendation

// POST /api/v1/ai/recommendations/similar-runs
getSimilarRuns(request: RecommendationRequest)
  → Array<SimilarRun>
```

#### B. 质量控制 (4个端点)

```typescript
// POST /api/v1/ai/qc/auto-analyze
runAutoQC(request: AutoQCRequest)
  → QCReport

// GET /api/v1/ai/qc/recommendations/{sampleId}
getQCRecommendations(sampleId: string)
  → string[]

// POST /api/v1/ai/qc/batch-analyze
batchQCAnalysis(sampleIds: string[])
  → QCReport[]

// POST /api/v1/ai/qc/predict-issues
predictQCIssues(metadata: Record<string, any>)
  → Array<{issue, probability, prevention}>
```

#### C. 异常检测 (4个端点)

```typescript
// POST /api/v1/ai/anomaly/detect
detectAnomalies(request: AnomalyDetectionRequest)
  → AnomalyDetectionReport

// GET /api/v1/ai/anomaly/{anomalyId}
getAnomalyDetails(anomalyId: string)
  → Anomaly

// POST /api/v1/ai/anomaly/{anomalyId}/fix
applyAnomalyFix(anomalyId: string, fixAction: string)
  → {success, message}

// WebSocket /ws/ai/anomaly/monitor/{projectId}
monitorAnomalies(projectId: string)
  → Real-time anomaly stream
```

#### D. 智能分组 (3个端点)

```typescript
// POST /api/v1/ai/grouping/smart-group
smartGroupSamples(request: SmartGroupingRequest)
  → SmartGroupingResult

// GET /api/v1/ai/grouping/suggest-comparisons/{projectId}
suggestComparisons(projectId: string)
  → Array<ComparisonSuggestion>

// POST /api/v1/ai/grouping/validate
validateGrouping(groups: Array<Group>)
  → {valid, issues, suggestions}
```

#### E. AI助手 (4个端点)

```typescript
// POST /api/v1/ai/assistant/conversations/{conversationId}/messages
sendAssistantMessage(conversationId, message, context)
  → AIAssistantMessage

// POST /api/v1/ai/assistant/conversations
createConversation(title, context)
  → AIAssistantConversation

// GET /api/v1/ai/assistant/conversations/{conversationId}
getConversation(conversationId)
  → AIAssistantConversation

// GET /api/v1/ai/assistant/conversations
listConversations()
  → AIAssistantConversation[]
```

#### F. 分析洞察 (3个端点)

```typescript
// GET /api/v1/ai/insights/project/{projectId}
getProjectInsights(projectId: string)
  → InsightsReport

// POST /api/v1/ai/insights/analyze
getAnalysisInsights(analysisType, data)
  → AnalysisInsight[]

// PUT /api/v1/ai/insights/{insightId}/reviewed
markInsightReviewed(insightId: string)
  → void
```

#### G. 预测功能 (3个端点)

```typescript
// POST /api/v1/ai/predictions/resources
predictResources(pipelineType, sampleCount, dataSize, parameters)
  → ResourcePrediction

// POST /api/v1/ai/predictions/success
predictSuccess(analysisConfig)
  → SuccessPrediction

// GET /api/v1/ai/predictions/timeline/{projectId}
predictTimeline(projectId: string)
  → Timeline
```

#### H. 系统状态 (2个端点)

```typescript
// GET /api/v1/ai/status
getSystemStatus()
  → AISystemStatus

// GET /api/v1/ai/available
isAvailable()
  → boolean
```

#### I. 反馈学习 (2个端点)

```typescript
// POST /api/v1/ai/feedback/{recommendationId}
submitFeedback(recommendationId, feedback)
  → void

// POST /api/v1/ai/feedback/prediction/{predictionId}
reportIncorrectPrediction(predictionId, actualOutcome, comments)
  → void
```

**总计**: 28个API端点需要实现

---

### 2. Analytics API (`analytics.service.ts`)

**方法数**: ~10个方法
**复杂度**: 中
**优先级**: P1 (核心功能)

```typescript
// GET /api/v1/analytics/overview
getOverview()
  → OverviewStats

// GET /api/v1/analytics/projects/{projectId}
getProjectAnalytics(projectId: string)
  → ProjectAnalytics

// GET /api/v1/analytics/samples/{sampleId}
getSampleAnalytics(sampleId: string)
  → SampleAnalytics

// GET /api/v1/analytics/tasks/statistics
getTaskStatistics(dateRange)
  → TaskStatistics

// GET /api/v1/analytics/pipelines/performance
getPipelinePerformance()
  → PipelinePerformance

// GET /api/v1/analytics/storage
getStorageAnalytics()
  → StorageAnalytics

// GET /api/v1/analytics/users/activity
getUserActivity(userId, dateRange)
  → UserActivity

// GET /api/v1/analytics/trends
getTrends(metric, period)
  → TrendData

// POST /api/v1/analytics/export
exportAnalytics(request: ExportRequest)
  → File
```

**总计**: 9个API端点需要实现

---

### 3. Notification API (`notification.service.ts`)

**方法数**: ~8个方法
**复杂度**: 低
**优先级**: P2 (重要功能)

```typescript
// GET /api/v1/notifications
getNotifications(params)
  → PaginatedResponse<Notification>

// GET /api/v1/notifications/{id}
getNotification(id: string)
  → Notification

// PUT /api/v1/notifications/{id}/read
markAsRead(id: string)
  → Notification

// PUT /api/v1/notifications/read-all
markAllAsRead()
  → {count}

// DELETE /api/v1/notifications/{id}
deleteNotification(id: string)
  → void

// GET /api/v1/notifications/unread/count
getUnreadCount()
  → {count}

// GET /api/v1/notifications/settings
getSettings()
  → NotificationSettings

// PUT /api/v1/notifications/settings
updateSettings(settings)
  → NotificationSettings

// WebSocket /ws/notifications
subscribeToNotifications()
  → Real-time notification stream
```

**总计**: 9个API端点需要实现

---

### 4. Admin API (`admin.service.ts` + `admin.enhanced.service.ts`)

**方法数**: ~15个方法
**复杂度**: 中
**优先级**: P2 (重要功能)

```typescript
// GET /api/v1/admin/dashboard
getDashboardStats()
  → AdminDashboardStats

// GET /api/v1/admin/users
getUsers(params)
  → PaginatedResponse<User>

// PUT /api/v1/admin/users/{id}
updateUser(id: string, data)
  → User

// DELETE /api/v1/admin/users/{id}
deleteUser(id: string)
  → void

// PUT /api/v1/admin/users/{id}/status
updateUserStatus(id: string, status)
  → User

// GET /api/v1/admin/system/health
getSystemHealth()
  → SystemHealth

// GET /api/v1/admin/system/logs
getLogs(params)
  → PaginatedResponse<LogEntry>

// GET /api/v1/admin/system/metrics
getSystemMetrics()
  → SystemMetrics

// POST /api/v1/admin/system/backup
createBackup()
  → BackupInfo

// GET /api/v1/admin/system/backups
listBackups()
  → Backup[]

// POST /api/v1/admin/system/restore
restoreBackup(backupId)
  → {success, message}

// GET /api/v1/admin/settings
getSettings()
  → SystemSettings

// PUT /api/v1/admin/settings
updateSettings(settings)
  → SystemSettings

// GET /api/v1/admin/audit-logs
getAuditLogs(params)
  → PaginatedResponse<AuditLog>
```

**总计**: 14个API端点需要实现

---

### 5. Stats API (`stats.service.ts`)

**方法数**: ~5个方法
**复杂度**: 低
**优先级**: P2 (重要功能)

```typescript
// GET /api/v1/stats/summary
getSummary()
  → StatsSummary

// GET /api/v1/stats/projects
getProjectStats()
  → ProjectStats

// GET /api/v1/stats/samples
getSampleStats()
  → SampleStats

// GET /api/v1/stats/tasks
getTaskStats()
  → TaskStats

// GET /api/v1/stats/storage
getStorageStats()
  → StorageStats
```

**总计**: 5个API端点需要实现

---

## 📋 实现优先级规划

### Phase 1: 核心功能API (P0 - 必须) ⚠️

**目标**: 确保基本CRUD功能正常工作

1. **验证现有API** (1天)
   - [ ] 测试auth API (login, register, token)
   - [ ] 测试users API (CRUD)
   - [ ] 测试projects API (CRUD)
   - [ ] 测试samples API (CRUD)
   - [ ] 测试files API (CRUD + upload)
   - [ ] 测试tasks API (CRUD + status)
   - [ ] 测试pipelines API (CRUD + templates)
   - [ ] 测试results API (CRUD)
   - [ ] 测试websocket API (connect)

2. **修复现有API问题** (1-2天)
   - [ ] 修复发现的bugs
   - [ ] 完善错误处理
   - [ ] 添加缺失的端点
   - [ ] 优化性能

### Phase 2: 统计和通知API (P1 - 高优先级) 📊

**目标**: 实现Dashboard和通知功能

3. **Stats API** (1天)
   - [ ] 实现统计数据聚合
   - [ ] 实现缓存机制
   - [ ] 创建统计service
   - [ ] 创建stats router
   - [ ] 测试所有端点

4. **Notification API** (1天)
   - [ ] 创建notification model
   - [ ] 实现通知CRUD
   - [ ] 实现WebSocket推送
   - [ ] 实现通知设置
   - [ ] 测试实时推送

### Phase 3: 分析功能API (P1 - 高优先级) 📈

**目标**: 实现数据分析和可视化

5. **Analytics API** (2天)
   - [ ] 实现数据聚合查询
   - [ ] 实现趋势分析
   - [ ] 实现导出功能
   - [ ] 创建analytics service
   - [ ] 创建analytics router
   - [ ] 测试所有分析端点

### Phase 4: 管理功能API (P2 - 中优先级) 👨‍💼

**目标**: 实现管理后台功能

6. **Admin API** (2天)
   - [ ] 实现用户管理端点
   - [ ] 实现系统健康检查
   - [ ] 实现日志查询
   - [ ] 实现系统设置
   - [ ] 实现审计日志
   - [ ] 创建admin service
   - [ ] 创建admin router
   - [ ] 添加权限控制

### Phase 5: AI功能API (P3 - 后续实现) 🤖

**目标**: 实现AI智能功能（需要AI模型支持）

7. **AI API - 基础框架** (3天)
   - [ ] 创建AI模型基础架构
   - [ ] 实现参数推荐框架
   - [ ] 实现QC检测框架
   - [ ] 实现异常检测框架
   - [ ] 创建ai service
   - [ ] 创建ai router
   - [ ] 返回模拟数据（待模型训练）

8. **AI API - 高级功能** (持续)
   - [ ] 集成真实AI模型
   - [ ] 训练参数推荐模型
   - [ ] 训练QC检测模型
   - [ ] 训练异常检测模型
   - [ ] 训练智能分组模型
   - [ ] 实现AI助手（可选LLM集成）
   - [ ] 优化预测性能

---

## 🔧 技术实现方案

### 1. 项目结构

```
backend/app/
├── api/v1/
│   ├── auth.py ✅
│   ├── users.py ✅
│   ├── projects.py ✅
│   ├── samples.py ✅
│   ├── files.py ✅
│   ├── tasks.py ✅
│   ├── pipelines.py ✅
│   ├── results.py ✅
│   ├── websocket.py ✅
│   ├── stats.py ❌ (待创建)
│   ├── notifications.py ❌ (待创建)
│   ├── analytics.py ❌ (待创建)
│   ├── admin.py ❌ (待创建)
│   └── ai.py ❌ (待创建)
├── models/
│   ├── notification.py ❌ (待创建)
│   ├── audit_log.py ❌ (待创建)
│   └── ai_model.py ❌ (待创建)
├── schemas/
│   ├── notification.py ❌ (待创建)
│   ├── analytics.py ❌ (待创建)
│   ├── admin.py ❌ (待创建)
│   └── ai.py ❌ (待创建)
├── services/
│   ├── stats_service.py ❌ (待创建)
│   ├── notification_service.py ❌ (待创建)
│   ├── analytics_service.py ❌ (待创建)
│   ├── admin_service.py ❌ (待创建)
│   └── ai_service.py ❌ (待创建)
└── utils/
    ├── ai_helpers.py ❌ (待创建)
    └── analytics_helpers.py ❌ (待创建)
```

### 2. 数据库Model设计

#### Notification Model

```python
class Notification(Base):
    __tablename__ = "notifications"

    id = Column(UUID, primary_key=True, default=uuid4)
    user_id = Column(UUID, ForeignKey("users.id"), nullable=False)
    type = Column(String, nullable=False)  # info, warning, error, success
    title = Column(String, nullable=False)
    message = Column(Text, nullable=False)
    data = Column(JSON, nullable=True)
    read = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.utcnow)

    user = relationship("User", back_populates="notifications")
```

#### AuditLog Model

```python
class AuditLog(Base):
    __tablename__ = "audit_logs"

    id = Column(UUID, primary_key=True, default=uuid4)
    user_id = Column(UUID, ForeignKey("users.id"), nullable=True)
    action = Column(String, nullable=False)
    resource_type = Column(String, nullable=False)
    resource_id = Column(String, nullable=True)
    details = Column(JSON, nullable=True)
    ip_address = Column(String, nullable=True)
    user_agent = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)

    user = relationship("User")
```

#### AIModel (可选，用于存储AI模型元数据)

```python
class AIModel(Base):
    __tablename__ = "ai_models"

    id = Column(UUID, primary_key=True, default=uuid4)
    name = Column(String, nullable=False)
    type = Column(String, nullable=False)  # recommendation, qc, anomaly, grouping
    version = Column(String, nullable=False)
    model_path = Column(String, nullable=False)
    accuracy = Column(Float, nullable=True)
    parameters = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, onupdate=datetime.utcnow)
```

### 3. Service层设计

所有新增的service应该遵循现有模式：

```python
# services/stats_service.py
class StatsService:
    def __init__(self, db: Session):
        self.db = db

    def get_summary(self) -> StatsSummary:
        """Get overall statistics summary"""
        pass

    def get_project_stats(self) -> ProjectStats:
        """Get project statistics"""
        pass

    # ... 其他方法
```

### 4. Router设计

所有新增的router应该遵循现有模式：

```python
# api/v1/stats.py
from fastapi import APIRouter, Depends
from app.core.deps import get_db, get_current_user
from app.services.stats_service import StatsService

router = APIRouter()

@router.get("/summary")
async def get_summary(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """Get statistics summary"""
    service = StatsService(db)
    return service.get_summary()
```

---

## 🧪 测试策略

### 单元测试

为每个新service编写单元测试：
```python
# tests/test_stats_service.py
def test_get_summary(db_session, test_user):
    service = StatsService(db_session)
    summary = service.get_summary()
    assert summary is not None
    assert summary.total_projects >= 0
```

### API集成测试

为每个新endpoint编写集成测试：
```python
# tests/test_stats_api.py
def test_get_summary_endpoint(client, auth_headers):
    response = client.get("/api/v1/stats/summary", headers=auth_headers)
    assert response.status_code == 200
    assert "total_projects" in response.json()
```

### 前后端联调测试

1. 启动后端服务
2. 启动前端开发服务器
3. 测试每个API端点
4. 验证数据格式和响应
5. 测试错误处理

---

## 📊 实现时间估算

| Phase | 内容 | 预计时间 | 优先级 |
|-------|------|---------|--------|
| Phase 1 | 验证和修复现有API | 2-3天 | P0 |
| Phase 2 | Stats + Notification API | 2天 | P1 |
| Phase 3 | Analytics API | 2天 | P1 |
| Phase 4 | Admin API | 2天 | P2 |
| Phase 5 | AI API基础框架 | 3天 | P3 |
| Phase 5 | AI模型集成和优化 | 持续 | P3 |

**总计**: 核心功能约11-12天，AI功能持续开发

---

## 🚀 立即开始

### 第一步：验证现有API

```bash
# 启动后端服务
cd backend
pip install -r requirements.txt
python init_db.py  # 初始化数据库
python create_admin.py  # 创建管理员账号
uvicorn app.main:app --reload

# 访问API文档
open http://localhost:8000/api/v1/docs
```

### 第二步：测试现有端点

使用Postman或curl测试所有现有端点，记录问题。

### 第三步：开始实现新API

按照优先级顺序，从Phase 2开始实现。

---

## 📝 下一步行动

**立即可做**:
1. ✅ 启动后端服务并测试现有API
2. ✅ 创建Stats API (高优先级)
3. ✅ 创建Notification API (高优先级)
4. ✅ 前后端联调测试

**短期目标** (1-2周):
- 完成Phase 1-3（Stats, Notification, Analytics）
- 实现Dashboard所需的所有数据端点
- 完成基础的前后端联调

**中期目标** (2-4周):
- 完成Phase 4 (Admin API)
- 实现管理后台所有功能
- 完成完整的前后端集成测试

**长期目标** (持续):
- 实现AI功能框架
- 训练和集成AI模型
- 持续优化性能

---

**文档创建时间**: 2024-01-01
**维护者**: NGSmodule Backend Team
**状态**: 规划中
