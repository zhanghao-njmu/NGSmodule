# NGSmodule 全面重构与Web平台开发计划

## 📋 项目概述

NGSmodule是一个用于NGS测序数据上游处理的分析平台，目标是构建一个现代化、用户友好的生物信息学工作站，支持终端和Web UI两种访问方式。

**当前状态**: 22,702行代码，功能完整但存在50+个需要修复的问题
**目标用户**: 科研人员、学生、教师（非编程背景）
**核心价值**: 全自动化、AI辅助、现代化UI、完整的数据分析流程

---

## 🎯 项目目标

### 终端版本
- ✅ 保持现有批量处理能力
- 🔧 修复已识别的bug和问题
- 📈 提升性能和稳定性
- 📚 完善文档和测试

### Web平台版本
- 🌐 构建现代化Web UI
- 👥 支持多用户系统（用户/管理员）
- 🤖 集成AI辅助功能
- 📊 实时流程监控和可视化
- 🔒 数据安全和权限管理

---

## 🏗️ Web UI 架构设计方案

### 1. 技术栈选型

#### 前端技术栈
```
框架: React 18+ (或 Vue 3+)
状态管理: Redux Toolkit / Zustand
UI组件库: Ant Design / Material-UI / Chakra UI
数据可视化:
  - Plotly.js (交互式图表)
  - ECharts (复杂图表)
  - Cytoscape.js (网络图)
  - IGV.js (基因组浏览器)
构建工具: Vite / Next.js
样式: TailwindCSS + CSS Modules
类型检查: TypeScript
```

**前端目录结构**:
```
web-ui/
├── src/
│   ├── components/          # 通用组件
│   │   ├── Layout/         # 布局组件
│   │   ├── Charts/         # 图表组件
│   │   ├── DataTable/      # 数据表格
│   │   ├── FileUploader/   # 文件上传
│   │   └── Pipeline/       # 流程组件
│   ├── pages/              # 页面组件
│   │   ├── Dashboard/      # 仪表板
│   │   ├── Projects/       # 项目管理
│   │   ├── Analysis/       # 分析页面
│   │   ├── Results/        # 结果展示
│   │   └── Admin/          # 管理后台
│   ├── services/           # API服务
│   ├── hooks/              # 自定义Hooks
│   ├── utils/              # 工具函数
│   ├── store/              # 状态管理
│   └── types/              # TypeScript类型
```

#### 后端技术栈
```
框架: FastAPI (Python 3.10+) / Django REST Framework
任务队列: Celery + Redis
消息队列: RabbitMQ / Redis Streams
数据库:
  - PostgreSQL (主数据库)
  - MongoDB (日志和文档存储)
  - Redis (缓存和会话)
文件存储: MinIO / S3
认证授权: JWT + OAuth2
API文档: OpenAPI (Swagger)
监控: Prometheus + Grafana
```

**后端目录结构**:
```
backend/
├── app/
│   ├── api/                # API路由
│   │   ├── v1/
│   │   │   ├── auth.py    # 认证
│   │   │   ├── users.py   # 用户管理
│   │   │   ├── projects.py # 项目
│   │   │   ├── pipelines.py # 管道
│   │   │   ├── tasks.py   # 任务
│   │   │   └── admin.py   # 管理接口
│   ├── core/               # 核心功能
│   │   ├── config.py      # 配置
│   │   ├── security.py    # 安全
│   │   └── database.py    # 数据库
│   ├── models/             # 数据模型
│   ├── schemas/            # Pydantic模式
│   ├── services/           # 业务逻辑
│   ├── workers/            # Celery任务
│   └── utils/              # 工具函数
├── scripts/                # Shell脚本包装器
│   └── pipeline_wrapper.py
├── tests/                  # 测试
└── alembic/               # 数据库迁移
```

#### 容器化和部署
```
Docker + Docker Compose
Kubernetes (生产环境)
Nginx (反向代理)
Let's Encrypt (SSL证书)
```

---

### 2. 系统架构设计

```
┌─────────────────────────────────────────────────────────────┐
│                      前端 (React/Vue)                        │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐   │
│  │ 仪表板   │  │ 项目管理 │  │ 分析页面 │  │ 管理后台 │   │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘   │
└────────────────────┬────────────────────────────────────────┘
                     │ REST API / WebSocket
┌────────────────────┴────────────────────────────────────────┐
│                   API网关 (Nginx/Kong)                       │
│                   认证中间件 (JWT)                           │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────┴────────────────────────────────────────┐
│                  应用服务器 (FastAPI)                        │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐   │
│  │ 用户服务 │  │ 项目服务 │  │ 任务服务 │  │ 文件服务 │   │
│  └──────────┘  └──────────┘  └──────────┘  └──────────┘   │
└────────┬───────────┬──────────┬──────────────┬─────────────┘
         │           │          │              │
    ┌────┴───┐  ┌───┴────┐ ┌───┴────────┐ ┌──┴──────┐
    │ 数据库 │  │ Redis  │ │ 任务队列   │ │ 对象存储│
    │(Postgres│  │(缓存)  │ │(Celery)    │ │ (MinIO) │
    └────────┘  └────────┘ └───┬────────┘ └─────────┘
                                │
                       ┌────────┴────────┐
                       │  NGS Pipeline   │
                       │   (Shell/R/Py)  │
                       └─────────────────┘
```

---

### 3. 数据库设计

#### 核心表结构

```sql
-- 用户表
CREATE TABLE users (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    username VARCHAR(50) UNIQUE NOT NULL,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    full_name VARCHAR(100),
    role VARCHAR(20) DEFAULT 'user', -- 'user' or 'admin'
    organization VARCHAR(100),
    is_active BOOLEAN DEFAULT true,
    storage_quota BIGINT DEFAULT 107374182400, -- 100GB
    storage_used BIGINT DEFAULT 0,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 项目表
CREATE TABLE projects (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES users(id) ON DELETE CASCADE,
    name VARCHAR(100) NOT NULL,
    description TEXT,
    project_type VARCHAR(20), -- 'rna-seq', 'dna-seq', 'sc-rna-seq'
    status VARCHAR(20) DEFAULT 'active', -- 'active', 'archived', 'deleted'
    config JSONB, -- 项目配置
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(user_id, name)
);

-- 样本表
CREATE TABLE samples (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    project_id UUID REFERENCES projects(id) ON DELETE CASCADE,
    sample_id VARCHAR(100) NOT NULL,
    run_id VARCHAR(100),
    group_name VARCHAR(50),
    layout VARCHAR(10), -- 'PE' or 'SE'
    batch_id VARCHAR(50),
    metadata JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(project_id, sample_id)
);

-- 文件表
CREATE TABLE files (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    sample_id UUID REFERENCES samples(id) ON DELETE CASCADE,
    filename VARCHAR(255) NOT NULL,
    file_path TEXT NOT NULL,
    file_type VARCHAR(20), -- 'fastq', 'bam', 'vcf', etc.
    file_size BIGINT,
    md5_checksum VARCHAR(32),
    upload_status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'uploading', 'completed', 'failed'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 流程任务表
CREATE TABLE pipeline_tasks (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    project_id UUID REFERENCES projects(id) ON DELETE CASCADE,
    task_name VARCHAR(100) NOT NULL,
    task_type VARCHAR(50), -- 'preAlignmentQC', 'Alignment', etc.
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'running', 'completed', 'failed', 'cancelled'
    progress FLOAT DEFAULT 0.0, -- 0-100
    started_at TIMESTAMP,
    completed_at TIMESTAMP,
    error_message TEXT,
    config JSONB,
    celery_task_id VARCHAR(255),
    log_file_path TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 结果表
CREATE TABLE results (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    task_id UUID REFERENCES pipeline_tasks(id) ON DELETE CASCADE,
    result_type VARCHAR(50), -- 'qc_report', 'alignment', 'quantification', 'de_analysis'
    result_path TEXT,
    metadata JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 系统日志表
CREATE TABLE system_logs (
    id BIGSERIAL PRIMARY KEY,
    user_id UUID REFERENCES users(id),
    action VARCHAR(50), -- 'login', 'create_project', 'run_pipeline', etc.
    resource_type VARCHAR(50),
    resource_id UUID,
    details JSONB,
    ip_address INET,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- 任务队列表 (用于优先级调度)
CREATE TABLE task_queue (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES users(id),
    priority INTEGER DEFAULT 10, -- 1-100, 越小优先级越高
    queue_position INTEGER,
    estimated_duration INTEGER, -- 预估分钟数
    status VARCHAR(20) DEFAULT 'queued',
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### MongoDB集合设计 (日志和文档存储)

```javascript
// 详细任务日志
{
  _id: ObjectId,
  task_id: "uuid",
  timestamp: ISODate,
  level: "INFO|WARNING|ERROR",
  message: "string",
  context: {
    sample_id: "string",
    step: "string",
    command: "string"
  }
}

// 分析报告
{
  _id: ObjectId,
  project_id: "uuid",
  report_type: "qc|alignment|de_analysis",
  content: {}, // 报告内容
  generated_at: ISODate
}
```

---

### 4. API设计 (RESTful)

#### 认证和用户管理

```
POST   /api/v1/auth/register          # 用户注册
POST   /api/v1/auth/login             # 用户登录
POST   /api/v1/auth/logout            # 用户登出
POST   /api/v1/auth/refresh           # 刷新Token
GET    /api/v1/users/me               # 获取当前用户信息
PUT    /api/v1/users/me               # 更新用户信息
GET    /api/v1/users/{user_id}        # 获取用户信息(管理员)
```

#### 项目管理

```
GET    /api/v1/projects               # 获取项目列表
POST   /api/v1/projects               # 创建项目
GET    /api/v1/projects/{id}          # 获取项目详情
PUT    /api/v1/projects/{id}          # 更新项目
DELETE /api/v1/projects/{id}          # 删除项目
POST   /api/v1/projects/{id}/archive  # 归档项目
```

#### 样本管理

```
GET    /api/v1/projects/{id}/samples        # 获取样本列表
POST   /api/v1/projects/{id}/samples        # 添加样本
POST   /api/v1/projects/{id}/samples/batch  # 批量导入样本
PUT    /api/v1/samples/{id}                 # 更新样本信息
DELETE /api/v1/samples/{id}                 # 删除样本
```

#### 文件管理

```
POST   /api/v1/files/upload           # 上传文件(分块上传)
GET    /api/v1/files/{id}             # 获取文件信息
GET    /api/v1/files/{id}/download    # 下载文件
DELETE /api/v1/files/{id}             # 删除文件
POST   /api/v1/files/validate         # 验证文件完整性
```

#### 流程管理

```
POST   /api/v1/pipelines/run          # 运行流程
GET    /api/v1/pipelines/tasks        # 获取任务列表
GET    /api/v1/pipelines/tasks/{id}   # 获取任务详情
POST   /api/v1/pipelines/tasks/{id}/cancel  # 取消任务
GET    /api/v1/pipelines/tasks/{id}/logs    # 获取任务日志
GET    /api/v1/pipelines/tasks/{id}/status  # 获取任务状态
```

#### 结果和可视化

```
GET    /api/v1/results/{task_id}           # 获取分析结果
GET    /api/v1/results/{task_id}/plots     # 获取图表数据
GET    /api/v1/results/{task_id}/download  # 下载结果文件
POST   /api/v1/results/{task_id}/report    # 生成报告
```

#### 管理后台

```
GET    /api/v1/admin/users             # 管理员-用户列表
PUT    /api/v1/admin/users/{id}        # 管理员-编辑用户
GET    /api/v1/admin/stats             # 系统统计信息
GET    /api/v1/admin/logs              # 系统日志
GET    /api/v1/admin/resources         # 资源使用情况
POST   /api/v1/admin/maintenance       # 维护模式
```

#### WebSocket接口 (实时更新)

```
ws://api/v1/ws/tasks/{task_id}        # 任务进度更新
ws://api/v1/ws/system                 # 系统通知
```

---

### 5. 前端页面设计

#### 5.1 登录/注册页面
```
特点:
- 现代化设计，毛玻璃效果
- 支持邮箱/用户名登录
- OAuth2第三方登录(可选)
- 忘记密码功能
- 注册需要邮箱验证
```

#### 5.2 用户仪表板 (Dashboard)
```
布局:
┌────────────────────────────────────────────────────┐
│  [Logo]  Dashboard  Projects  Analysis   [User▼]  │
├────────────────────────────────────────────────────┤
│                                                     │
│  欢迎回来, [用户名]                                  │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ │
│  │ 项目总数 │ │ 运行中  │ │ 已完成  │ │ 存储使用│ │
│  │   12    │ │   3     │ │   45    │ │  65/100│ │
│  └─────────┘ └─────────┘ └─────────┘ └─────────┘ │
│                                                     │
│  最近项目                                [查看全部→] │
│  ┌─────────────────────────────────────────────┐  │
│  │ Project A    RNA-seq   运行中   70%  [详情]│  │
│  │ Project B    scRNA-seq 已完成   100% [详情]│  │
│  └─────────────────────────────────────────────┘  │
│                                                     │
│  运行中的任务                          [查看全部→]  │
│  ┌─────────────────────────────────────────────┐  │
│  │ Sample1-Alignment  ████████░░  80%  [日志]│  │
│  │ Sample2-QC         ████░░░░░  40%  [日志]│  │
│  └─────────────────────────────────────────────┘  │
└────────────────────────────────────────────────────┘
```

**组件**:
- 统计卡片 (项目数、任务状态、存储)
- 最近项目列表
- 运行中任务实时进度
- 快速操作按钮 (新建项目、上传数据)
- 系统通知

#### 5.3 项目管理页面
```
功能:
- 项目列表 (卡片/表格视图切换)
- 搜索和筛选 (按类型、状态、日期)
- 创建新项目向导
- 项目详情页
  - 样本管理
  - 配置编辑
  - 文件上传
  - 流程运行历史
```

**项目创建向导**:
```
步骤1: 选择分析类型
  ○ RNA-seq
  ○ DNA-seq (WGS/WES)
  ○ Single-cell RNA-seq
  ○ Bisulfite-seq
  ○ ATAC-seq

步骤2: 项目信息
  - 项目名称: [________]
  - 描述: [____________]
  - 物种: [Human▼]
  - 基因组版本: [GRCh38▼]

步骤3: 样本信息
  - CSV导入
  - 手动添加
  - 模板下载

步骤4: 流程配置
  - 预设模板选择
  - 自定义参数 (高级)
```

#### 5.4 数据上传页面
```
特点:
- 拖拽上传
- 分块上传 (大文件支持)
- 断点续传
- MD5校验
- 批量上传
- 进度可视化
- 上传历史
```

**界面**:
```
┌─────────────────────────────────────────────┐
│  拖拽文件到这里或点击选择                      │
│                                              │
│      ☁️                                      │
│   支持: .fastq, .fastq.gz, .fq.gz           │
│   单文件最大: 50GB                            │
│   批量上传: 无限制                            │
└─────────────────────────────────────────────┘

上传列表:
┌─────────────────────────────────────────────┐
│ sample1_R1.fastq.gz  ████████████  100%  ✓│
│ sample1_R2.fastq.gz  ████████░░░░   85%  ⟳│
│ sample2_R1.fastq.gz  ░░░░░░░░░░░░    0%  ⏸│
└─────────────────────────────────────────────┘
```

#### 5.5 流程运行页面
```
特点:
- 可视化流程图
- 实时进度监控
- 日志查看器
- 资源使用监控
- 错误诊断和建议
```

**流程图示例**:
```
Raw Data → QC → Trimming → Alignment → Quantification → DE Analysis
  ✓        ✓      ✓           ⟳             ○              ○
         80%
         预计剩余: 45分钟
```

**日志查看器**:
```
┌─────────────────────────────────────────────┐
│ [INFO] 2025-11-21 10:23:45                  │
│ Starting alignment for Sample1...            │
│                                              │
│ [INFO] Using HISAT2 aligner                 │
│ [INFO] Threads: 16                          │
│ [PROGRESS] Aligned: 15,234,567 reads (45%)  │
│                                              │
│ [实时更新]                     [下载日志]     │
└─────────────────────────────────────────────┘
```

#### 5.6 结果分析页面
```
功能:
- 多维度数据可视化
- 交互式图表 (Plotly)
- IGV基因组浏览器集成
- 结果下载
- 报告生成
- 结果分享
```

**可视化组件**:
- QC报告: FastQC结果、质量分布、GC含量
- 比对统计: 比对率、覆盖度、基因体分布
- 定量结果: 基因表达热图、PCA、相关性矩阵
- 差异分析: 火山图、MA图、基因集富集
- 单细胞: UMAP、tSNE、marker基因

#### 5.7 管理后台页面 (仅管理员)
```
功能:
- 用户管理 (创建、编辑、删除、权限)
- 系统监控 (CPU、内存、磁盘、任务队列)
- 配置管理 (全局配置、参考基因组)
- 任务调度 (优先级、资源分配)
- 日志查看 (系统日志、错误日志)
- 统计报表 (用户活跃度、资源使用)
```

**系统监控仪表板**:
```
┌─────────────────────────────────────────────┐
│  CPU使用率   [████████░░] 78%               │
│  内存使用    [██████░░░░] 65%               │
│  磁盘使用    [███░░░░░░░] 32%               │
│                                              │
│  当前运行任务: 12                            │
│  队列中任务: 8                               │
│  活跃用户: 45                                │
└─────────────────────────────────────────────┘
```

---

### 6. AI功能集成

#### 6.1 智能参数推荐
```python
# 基于历史成功任务推荐最优参数
def recommend_parameters(project_type, species, sample_info):
    """
    使用机器学习模型推荐参数
    - 数据: 历史项目的参数和成功率
    - 特征: 样本量、测序深度、数据质量
    - 输出: 推荐的参数配置和置信度
    """
    pass
```

#### 6.2 自动QC诊断
```python
def diagnose_qc_issues(qc_results):
    """
    自动识别数据质量问题
    - 低质量碱基百分比异常
    - GC含量偏差
    - 接头污染
    - rRNA污染
    输出: 问题诊断和修复建议
    """
    pass
```

#### 6.3 智能任务调度
```python
def optimize_task_schedule(pending_tasks, system_resources):
    """
    根据资源可用性和任务优先级优化调度
    - 考虑用户配额
    - 预测任务运行时间
    - 负载均衡
    """
    pass
```

#### 6.4 结果解读助手
```python
def interpret_results(analysis_type, results):
    """
    AI助手解读分析结果
    - 差异基因生物学意义
    - 通路富集解读
    - 文献推荐
    """
    pass
```

#### 6.5 自然语言查询
```
用户输入: "显示TP53基因在所有样本中的表达水平"
系统: 自动生成查询 → 可视化结果
```

---

### 7. 实时通信和通知系统

#### WebSocket实现

```python
# backend/app/websockets/task_monitor.py
from fastapi import WebSocket

async def task_progress_stream(websocket: WebSocket, task_id: str):
    """实时推送任务进度"""
    await websocket.accept()
    while True:
        progress = await get_task_progress(task_id)
        await websocket.send_json({
            "task_id": task_id,
            "progress": progress,
            "status": "running",
            "message": f"Processing {progress}%"
        })
        if progress >= 100:
            break
        await asyncio.sleep(2)
```

#### 通知类型
- 任务完成通知
- 任务失败告警
- 磁盘空间预警
- 系统维护通知
- 新功能上线通知

---

### 8. 安全和权限设计

#### 认证机制
```
- JWT Token (15分钟有效期)
- Refresh Token (7天有效期)
- 双因素认证 (可选)
- API Key (用于程序化访问)
```

#### 权限模型
```
角色:
- Admin: 全部权限
- User: 自己的项目和数据
- Guest: 只读权限 (公共项目)

资源级权限:
- Project: owner, collaborator, viewer
- Data: read, write, delete
- Pipeline: execute, cancel
```

#### 数据安全
```
- 传输加密: HTTPS/TLS
- 存储加密: 敏感数据加密存储
- 访问日志: 所有操作记录
- 数据隔离: 用户数据物理隔离
- 定期备份: 自动备份机制
```

---

## 🔧 终端版本重构计划

### Phase 1: 代码质量提升 (2-3周)

#### 1.1 添加全局错误处理
```bash
# lib/error_handler.sh
set -euo pipefail

trap 'error_handler $? $LINENO' ERR
trap 'cleanup; exit 130' INT TERM

error_handler() {
    local exit_code=$1
    local line_number=$2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Command failed with exit code $exit_code at line $line_number" >&2
    echo "Stack trace:" >&2
    local frame=0
    while caller $frame >&2; do
        ((frame++))
    done
    exit $exit_code
}

cleanup() {
    # 清理临时文件和进程
    echo "Cleaning up..." >&2
}
```

#### 1.2 参数化硬编码路径
```bash
# config/defaults.sh - 创建默认配置
: "${IGENOMES_DIR:=/data/reference/iGenomes}"
: "${FASTQ_SCREEN_CONFIG:=/data/reference/FastQ_Screen/fastq_screen.conf}"
: "${SORTMERNA_DIR:=/data/reference/SortmeRNA}"
: "${THREADS:=$(nproc)}"

export IGENOMES_DIR FASTQ_SCREEN_CONFIG SORTMERNA_DIR THREADS

# 在主脚本中
source "$(dirname "$0")/config/defaults.sh"
source "${CUSTOM_CONFIG:-$HOME/.ngsmodule/config.sh}"
```

#### 1.3 创建共享函数库
```bash
# lib/common.sh
require_file() {
    local file=$1
    [[ -f "$file" && -r "$file" ]] || {
        log_error "Required file not found or not readable: $file"
        return 1
    }
}

require_dir() {
    local dir=$1
    [[ -d "$dir" && -w "$dir" ]] || {
        log_error "Required directory not found or not writable: $dir"
        return 1
    }
}

require_command() {
    local cmd=$1
    command -v "$cmd" &>/dev/null || {
        log_error "Required command not found: $cmd"
        return 1
    }
}

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $*" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $*" >&2
}
```

#### 1.4 统一R包管理
```r
# lib/check_packages.R
check_and_install_packages <- function() {
    required_packages <- list(
        "dplyr" = "1.0.0",
        "limma" = "3.50.0",
        "edgeR" = "3.36.0",
        "DESeq2" = "1.34.0",
        "Seurat" = "4.3.0"
        # ... 所有依赖
    )

    missing_packages <- character()

    for (pkg in names(required_packages)) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing_packages <- c(missing_packages, pkg)
        } else {
            # 检查版本
            installed_version <- packageVersion(pkg)
            required_version <- required_packages[[pkg]]
            if (installed_version < required_version) {
                warning(sprintf(
                    "Package %s version %s is installed, but %s is required",
                    pkg, installed_version, required_version
                ))
            }
        }
    }

    if (length(missing_packages) > 0) {
        stop(sprintf(
            "Missing required packages: %s\nPlease install them using BiocManager::install()",
            paste(missing_packages, collapse = ", ")
        ))
    }

    invisible(TRUE)
}
```

### Phase 2: 功能增强 (4-6周)

#### 2.1 添加断点续传功能
```bash
# lib/checkpoint.sh
save_checkpoint() {
    local sample=$1
    local step=$2
    local status=$3
    echo "${sample},${step},${status},$(date +%s)" >> "$maindir/.checkpoints"
}

load_checkpoint() {
    local sample=$1
    local step=$2
    if [[ -f "$maindir/.checkpoints" ]]; then
        grep "^${sample},${step},completed," "$maindir/.checkpoints" &>/dev/null
        return $?
    fi
    return 1
}

should_skip_step() {
    local sample=$1
    local step=$2
    local force=${3:-FALSE}

    if [[ "$force" == "TRUE" ]]; then
        return 1  # 不跳过
    fi

    if load_checkpoint "$sample" "$step"; then
        log_info "Skipping $step for $sample (already completed)"
        return 0  # 跳过
    fi

    return 1  # 不跳过
}
```

#### 2.2 增强日志系统
```bash
# lib/logging.sh
LOG_LEVEL=${LOG_LEVEL:-INFO}  # DEBUG, INFO, WARNING, ERROR
LOG_FILE=${LOG_FILE:-/dev/stderr}

log() {
    local level=$1
    shift
    local message="$*"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

    # 日志级别数值
    declare -A log_levels=([DEBUG]=0 [INFO]=1 [WARNING]=2 [ERROR]=3)

    if (( ${log_levels[$level]} >= ${log_levels[$LOG_LEVEL]} )); then
        echo "[$timestamp] [$level] $message" >> "$LOG_FILE"
    fi
}

log_command() {
    local cmd="$*"
    log INFO "Executing: $cmd"
    eval "$cmd" 2>&1 | tee -a "$LOG_FILE"
    local exit_code=${PIPESTATUS[0]}
    if [[ $exit_code -ne 0 ]]; then
        log ERROR "Command failed with exit code $exit_code: $cmd"
        return $exit_code
    fi
    log INFO "Command completed successfully"
    return 0
}
```

#### 2.3 资源监控
```bash
# lib/resource_monitor.sh
monitor_resources() {
    local interval=${1:-60}  # 监控间隔(秒)
    local log_file=${2:-resource_usage.log}

    while true; do
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        local cpu_usage=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)
        local mem_usage=$(free | grep Mem | awk '{printf("%.2f", $3/$2 * 100.0)}')
        local disk_usage=$(df -h "$maindir" | tail -1 | awk '{print $5}' | cut -d'%' -f1)

        echo "$timestamp,$cpu_usage,$mem_usage,$disk_usage" >> "$log_file"
        sleep "$interval"
    done
}
```

### Phase 3: 测试和文档 (3-4周)

#### 3.1 创建测试框架
```bash
# tests/test_framework.sh
source "$(dirname "$0")/../lib/common.sh"

test_count=0
test_passed=0
test_failed=0

assert_equals() {
    local expected=$1
    local actual=$2
    local message=${3:-""}

    ((test_count++))

    if [[ "$expected" == "$actual" ]]; then
        ((test_passed++))
        echo "✓ Test $test_count passed: $message"
    else
        ((test_failed++))
        echo "✗ Test $test_count failed: $message"
        echo "  Expected: $expected"
        echo "  Actual: $actual"
    fi
}

assert_file_exists() {
    local file=$1
    local message=${2:-""}

    ((test_count++))

    if [[ -f "$file" ]]; then
        ((test_passed++))
        echo "✓ Test $test_count passed: $message"
    else
        ((test_failed++))
        echo "✗ Test $test_count failed: $message"
        echo "  File does not exist: $file"
    fi
}

run_tests() {
    echo "Running tests..."
    echo "================"

    # 运行所有测试
    for test_file in tests/test_*.sh; do
        echo "Running $test_file..."
        source "$test_file"
    done

    echo "================"
    echo "Tests: $test_count"
    echo "Passed: $test_passed"
    echo "Failed: $test_failed"

    if (( test_failed > 0 )); then
        exit 1
    fi
}
```

#### 3.2 示例测试
```bash
# tests/test_common.sh
source "$(dirname "$0")/test_framework.sh"
source "$(dirname "$0")/../lib/common.sh"

# 测试 require_file
test_require_file() {
    # 创建测试文件
    local test_file=$(mktemp)
    assert_equals 0 "$(require_file "$test_file"; echo $?)" "require_file should succeed for existing file"
    rm "$test_file"

    # 测试不存在的文件
    assert_equals 1 "$(require_file "/nonexistent/file"; echo $?)" "require_file should fail for missing file"
}

# 运行测试
test_require_file
```

#### 3.3 文档结构
```
docs/
├── README.md                  # 项目概述
├── installation.md            # 安装指南
├── quickstart.md              # 快速开始
├── user-guide/                # 用户指南
│   ├── rna-seq.md
│   ├── dna-seq.md
│   ├── single-cell.md
│   └── configuration.md
├── api/                       # API文档
│   ├── rest-api.md
│   └── websocket-api.md
├── developer-guide/           # 开发者指南
│   ├── architecture.md
│   ├── contributing.md
│   └── testing.md
└── troubleshooting.md         # 故障排除
```

---

## 📦 部署方案

### Docker Compose部署 (开发/小规模)

```yaml
# docker-compose.yml
version: '3.8'

services:
  # 前端
  frontend:
    build: ./web-ui
    ports:
      - "3000:3000"
    environment:
      - REACT_APP_API_URL=http://localhost:8000
    depends_on:
      - backend

  # 后端API
  backend:
    build: ./backend
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql://user:pass@postgres:5432/ngsmodule
      - REDIS_URL=redis://redis:6379/0
      - CELERY_BROKER_URL=redis://redis:6379/0
    volumes:
      - ./data:/data
      - ./NGSmodule:/app/NGSmodule
    depends_on:
      - postgres
      - redis
      - minio

  # Celery工作节点
  celery-worker:
    build: ./backend
    command: celery -A app.workers worker --loglevel=info --concurrency=4
    environment:
      - DATABASE_URL=postgresql://user:pass@postgres:5432/ngsmodule
      - REDIS_URL=redis://redis:6379/0
    volumes:
      - ./data:/data
      - ./NGSmodule:/app/NGSmodule
    depends_on:
      - postgres
      - redis

  # Celery监控
  flower:
    build: ./backend
    command: celery -A app.workers flower
    ports:
      - "5555:5555"
    environment:
      - CELERY_BROKER_URL=redis://redis:6379/0
    depends_on:
      - redis
      - celery-worker

  # PostgreSQL数据库
  postgres:
    image: postgres:15
    environment:
      - POSTGRES_USER=user
      - POSTGRES_PASSWORD=pass
      - POSTGRES_DB=ngsmodule
    volumes:
      - postgres_data:/var/lib/postgresql/data
    ports:
      - "5432:5432"

  # Redis缓存
  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"

  # MinIO对象存储
  minio:
    image: minio/minio
    command: server /data --console-address ":9001"
    environment:
      - MINIO_ROOT_USER=minioadmin
      - MINIO_ROOT_PASSWORD=minioadmin
    ports:
      - "9000:9000"
      - "9001:9001"
    volumes:
      - minio_data:/data

  # Nginx反向代理
  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
      - ./ssl:/etc/nginx/ssl
    depends_on:
      - frontend
      - backend

volumes:
  postgres_data:
  minio_data:
```

### Kubernetes部署 (生产环境)

```yaml
# k8s/deployment.yaml
apiVersion: apps/v1
kind: Deployment
metadata:
  name: ngsmodule-backend
spec:
  replicas: 3
  selector:
    matchLabels:
      app: ngsmodule-backend
  template:
    metadata:
      labels:
        app: ngsmodule-backend
    spec:
      containers:
      - name: backend
        image: ngsmodule/backend:latest
        ports:
        - containerPort: 8000
        env:
        - name: DATABASE_URL
          valueFrom:
            secretKeyRef:
              name: ngsmodule-secrets
              key: database-url
        resources:
          requests:
            memory: "2Gi"
            cpu: "1000m"
          limits:
            memory: "4Gi"
            cpu: "2000m"
        livenessProbe:
          httpGet:
            path: /health
            port: 8000
          initialDelaySeconds: 30
          periodSeconds: 10
        readinessProbe:
          httpGet:
            path: /ready
            port: 8000
          initialDelaySeconds: 5
          periodSeconds: 5
---
apiVersion: v1
kind: Service
metadata:
  name: ngsmodule-backend
spec:
  selector:
    app: ngsmodule-backend
  ports:
  - port: 8000
    targetPort: 8000
  type: ClusterIP
```

---

## 📊 项目时间线

### 总体规划 (6-9个月)

```
阶段1: 基础设施和架构 (Week 1-4)
├─ Week 1-2: 后端框架搭建
│  ├─ FastAPI项目初始化
│  ├─ 数据库设计和迁移
│  ├─ 认证授权系统
│  └─ 基础API实现
└─ Week 3-4: 前端框架搭建
   ├─ React项目初始化
   ├─ 路由和状态管理
   ├─ UI组件库集成
   └─ 登录/注册页面

阶段2: 核心功能开发 (Week 5-16)
├─ Week 5-7: 项目和样本管理
│  ├─ 项目CRUD API
│  ├─ 样本管理API
│  ├─ 前端项目管理页面
│  └─ 样本信息导入
├─ Week 8-10: 文件上传系统
│  ├─ 分块上传API
│  ├─ MinIO集成
│  ├─ 前端上传组件
│  └─ 进度监控
├─ Week 11-13: 流程集成
│  ├─ Shell脚本包装器
│  ├─ Celery任务调度
│  ├─ 流程执行API
│  └─ WebSocket实时更新
└─ Week 14-16: 结果展示
   ├─ 结果存储和检索
   ├─ 可视化组件
   ├─ 交互式图表
   └─ 报告生成

阶段3: 终端版本重构 (Week 5-12，并行)
├─ Week 5-7: 代码质量提升
│  ├─ 错误处理增强
│  ├─ 参数化配置
│  ├─ 共享函数库
│  └─ 日志系统改进
├─ Week 8-10: 功能增强
│  ├─ 断点续传
│  ├─ 资源监控
│  ├─ 性能优化
│  └─ 依赖管理
└─ Week 11-12: 测试和文档
   ├─ 测试框架
   ├─ 单元测试
   ├─ 集成测试
   └─ 文档编写

阶段4: 高级功能 (Week 17-24)
├─ Week 17-19: AI功能
│  ├─ 参数推荐模型
│  ├─ QC诊断系统
│  ├─ 智能调度
│  └─ 结果解读
├─ Week 20-21: 管理后台
│  ├─ 用户管理
│  ├─ 系统监控
│  ├─ 资源管理
│  └─ 统计报表
└─ Week 22-24: 协作功能
   ├─ 项目分享
   ├─ 结果导出
   ├─ 评论和注释
   └─ 版本控制

阶段5: 测试和优化 (Week 25-30)
├─ Week 25-27: 全面测试
│  ├─ 功能测试
│  ├─ 性能测试
│  ├─ 安全测试
│  └─ 用户测试
├─ Week 28-29: 优化和修复
│  ├─ 性能优化
│  ├─ Bug修复
│  ├─ UI/UX改进
│  └─ 文档完善
└─ Week 30: 发布准备
   ├─ 部署文档
   ├─ 用户手册
   ├─ 培训材料
   └─ 上线检查

阶段6: 发布和迭代 (Week 31+)
├─ Beta测试
├─ 收集反馈
├─ 快速迭代
└─ 正式发布
```

---

## 🎨 UI/UX设计原则

### 设计理念
1. **简洁优雅**: 现代化、扁平化设计
2. **直观易用**: 符合用户直觉，学习成本低
3. **视觉层次**: 清晰的信息架构
4. **响应式**: 支持桌面、平板、移动端
5. **可访问性**: 符合WCAG 2.1标准

### 色彩方案
```css
/* 主题色 */
--primary-color: #2196F3;      /* 蓝色 - 科技感 */
--secondary-color: #00BCD4;     /* 青色 - 生物感 */
--accent-color: #FF4081;        /* 粉色 - 强调 */
--success-color: #4CAF50;       /* 绿色 */
--warning-color: #FF9800;       /* 橙色 */
--error-color: #F44336;         /* 红色 */
--info-color: #2196F3;          /* 蓝色 */

/* 中性色 */
--background-light: #FAFAFA;
--background-dark: #121212;
--text-primary: #212121;
--text-secondary: #757575;
--border-color: #E0E0E0;

/* 数据可视化色板 */
--data-colors: ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"];
```

### 字体
```css
/* 英文 */
font-family: 'Inter', 'Roboto', system-ui, sans-serif;

/* 中文 */
font-family: 'Noto Sans SC', 'PingFang SC', 'Microsoft YaHei', sans-serif;

/* 代码 */
font-family: 'Fira Code', 'Consolas', 'Monaco', monospace;
```

### 动画和交互
```css
/* 过渡动画 */
transition: all 0.3s cubic-bezier(0.4, 0.0, 0.2, 1);

/* 悬停效果 */
.card:hover {
    transform: translateY(-4px);
    box-shadow: 0 8px 16px rgba(0,0,0,0.1);
}

/* 加载动画 */
@keyframes pulse {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.5; }
}
```

---

## 📈 成功指标

### 技术指标
- **页面加载时间**: < 2秒
- **API响应时间**: < 200ms (P95)
- **系统可用性**: > 99.5%
- **数据处理吞吐量**: > 100 样本/天
- **并发用户**: > 100

### 用户指标
- **注册用户**: > 1000 (第一年)
- **活跃用户**: > 500 (月活)
- **用户留存率**: > 60% (月留存)
- **任务成功率**: > 95%
- **用户满意度**: > 4.5/5

### 业务指标
- **项目数量**: > 5000 (第一年)
- **数据处理量**: > 10TB
- **平均任务时长**: 减少30%
- **错误率**: < 5%

---

## 🔄 持续改进计划

### 短期 (3个月)
- 收集用户反馈
- 修复高优先级Bug
- 性能优化
- 文档完善

### 中期 (6个月)
- 新增分析类型 (ChIP-seq, ATAC-seq)
- 增强AI功能
- 移动端App
- 国际化支持

### 长期 (12个月)
- 云服务集成 (AWS, Azure, 阿里云)
- 机器学习模型训练平台
- 协作分析功能
- 知识库和社区

---

## 📚 参考资源

### 类似平台
- Galaxy Project
- Terra/Broad
- DNAnexus
- BaseSpace (Illumina)

### 技术文档
- FastAPI: https://fastapi.tiangolo.com/
- React: https://react.dev/
- Plotly: https://plotly.com/javascript/
- Celery: https://docs.celeryq.dev/

### 生物信息学资源
- Bioconductor: https://bioconductor.org/
- Biostars: https://www.biostars.org/
- SEQanswers: http://seqanswers.com/

---

## 🤝 团队协作

### 推荐团队配置
```
项目经理 (PM): 1人
前端工程师: 2-3人
后端工程师: 2-3人
DevOps工程师: 1人
生物信息学专家: 1-2人
UI/UX设计师: 1人
测试工程师: 1人
技术文档工程师: 1人
```

### 开发流程
```
1. Sprint计划 (每2周)
2. 每日站会 (15分钟)
3. 代码审查 (所有PR)
4. 集成测试 (每天)
5. Sprint回顾 (每2周)
```

### Git工作流
```
main (生产分支)
  ↑
develop (开发分支)
  ↑
feature/* (功能分支)
hotfix/* (热修复分支)
```

---

## 📄 许可证和开源

建议采用 **MIT License** 或 **Apache License 2.0**，允许：
- 商业使用
- 修改
- 分发
- 私有使用

---

## 联系和支持

- 项目仓库: https://github.com/[org]/NGSmodule
- 文档: https://docs.ngsmodule.org
- 社区: https://community.ngsmodule.org
- 邮箱: support@ngsmodule.org

---

**最后更新**: 2025-11-21
**版本**: 1.0.0
