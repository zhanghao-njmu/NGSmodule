# NGSmodule 部署指南

## 🚀 快速开始

### 前提条件

- Docker >= 20.10
- Docker Compose >= 2.0
- 至少 4GB 可用内存
- 至少 10GB 可用磁盘空间

### 一键启动

```bash
# 克隆项目（如果还没有）
git clone <repository-url>
cd NGSmodule

# 运行快速启动脚本
./start.sh
```

启动脚本会自动：
1. 检查 Docker 和 Docker Compose 是否安装
2. 创建环境变量文件（从 .env.example 复制）
3. 启动所有服务（PostgreSQL, Redis, MinIO, Backend, Celery, Frontend）
4. 初始化数据库表
5. 创建默认管理员账户

### 手动启动

如果你想手动控制每个步骤：

#### 1. 配置环境变量

```bash
# 后端环境变量
cp backend/.env.example backend/.env

# 前端环境变量
cp frontend/.env.example frontend/.env

# 根据需要修改配置
vim backend/.env
vim frontend/.env
```

#### 2. 启动服务

```bash
# 启动所有服务
docker-compose up -d

# 查看服务状态
docker-compose ps

# 查看日志
docker-compose logs -f
```

#### 3. 初始化数据库

```bash
# 创建数据库表
docker-compose exec backend python init_db.py

# 创建管理员账户
docker-compose exec backend python create_admin.py
```

---

## 📍 访问应用

启动成功后，可以通过以下地址访问：

| 服务 | 地址 | 说明 |
|------|------|------|
| **前端应用** | http://localhost:3000 | 主Web界面 |
| **后端API** | http://localhost:8000 | RESTful API |
| **API文档** | http://localhost:8000/api/v1/docs | Swagger UI |
| **Flower** | http://localhost:5555 | Celery任务监控 |
| **MinIO控制台** | http://localhost:9001 | 对象存储管理 |

### 默认登录凭据

```
用户名: admin
密码: admin123
```

⚠️ **重要**: 首次登录后请立即修改密码！

---

## 🗄️ 服务架构

```
┌─────────────────────────────────────────────┐
│          Frontend (React)                    │
│          Port: 3000                          │
└──────────────────┬──────────────────────────┘
                   │ HTTP/WebSocket
┌──────────────────▼──────────────────────────┐
│          Backend (FastAPI)                   │
│          Port: 8000                          │
└───┬───────────┬──────────────┬──────────────┘
    │           │              │
    ↓           ↓              ↓
┌────────┐  ┌────────┐  ┌──────────┐
│PostgreSQL│ │ Redis  │ │  MinIO   │
│Port:5432│  │Port:6379│ │Port:9000│
└────────┘  └────┬───┘  └──────────┘
                 │
            ┌────▼─────┐
            │  Celery  │
            │  Worker  │
            └──────────┘
```

---

## 🛠️ 开发模式

### 本地开发（无 Docker）

#### 后端开发

```bash
cd backend

# 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# 安装依赖
pip install -r requirements.txt

# 配置环境变量
cp .env.example .env
vim .env  # 修改DATABASE_URL等配置

# 启动数据库服务（使用Docker）
docker-compose up -d postgres redis minio

# 初始化数据库
python init_db.py
python create_admin.py

# 启动后端服务
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000

# 另开终端，启动Celery Worker
celery -A app.workers.celery_app worker --loglevel=info
```

#### 前端开发

```bash
cd frontend

# 安装依赖
npm install

# 配置环境变量
cp .env.example .env

# 启动开发服务器
npm run dev
```

---

## 🔧 常用命令

### Docker Compose 命令

```bash
# 启动所有服务
docker-compose up -d

# 停止所有服务
docker-compose down

# 停止并删除所有数据（包括数据库）
docker-compose down -v

# 查看服务状态
docker-compose ps

# 查看日志
docker-compose logs -f

# 查看特定服务日志
docker-compose logs -f backend
docker-compose logs -f frontend

# 重启服务
docker-compose restart backend

# 进入容器
docker-compose exec backend bash
docker-compose exec frontend sh

# 重新构建镜像
docker-compose build --no-cache
```

### 数据库管理

```bash
# 连接到PostgreSQL
docker-compose exec postgres psql -U ngsmodule -d ngsmodule

# 备份数据库
docker-compose exec postgres pg_dump -U ngsmodule ngsmodule > backup.sql

# 恢复数据库
docker-compose exec -T postgres psql -U ngsmodule ngsmodule < backup.sql
```

### 后端命令

```bash
# 创建数据库迁移
docker-compose exec backend alembic revision --autogenerate -m "description"

# 执行数据库迁移
docker-compose exec backend alembic upgrade head

# 创建管理员用户
docker-compose exec backend python create_admin.py

# 进入Python Shell
docker-compose exec backend python
```

---

## 📦 数据持久化

项目使用 Docker Volumes 持久化数据：

| Volume | 用途 |
|--------|------|
| `postgres_data` | PostgreSQL 数据库文件 |
| `redis_data` | Redis 持久化数据 |
| `minio_data` | MinIO 对象存储文件 |
| `upload_data` | 用户上传的文件 |
| `work_data` | NGS 分析工作目录 |
| `result_data` | NGS 分析结果 |

### 备份数据

```bash
# 导出所有数据卷
docker run --rm \
  -v ngsmodule_postgres_data:/data/postgres \
  -v ngsmodule_minio_data:/data/minio \
  -v $(pwd):/backup \
  alpine tar czf /backup/ngsmodule-backup-$(date +%Y%m%d).tar.gz /data
```

---

## 🔍 故障排除

### 服务无法启动

```bash
# 查看详细日志
docker-compose logs -f

# 检查端口是否被占用
lsof -i :3000  # 前端
lsof -i :8000  # 后端
lsof -i :5432  # PostgreSQL

# 清理并重新启动
docker-compose down -v
docker-compose up -d
```

### 数据库连接失败

```bash
# 检查PostgreSQL是否运行
docker-compose ps postgres

# 测试数据库连接
docker-compose exec backend python -c "from app.core.database import engine; print(engine.connect())"

# 重启PostgreSQL
docker-compose restart postgres
```

### 前端无法连接后端

1. 检查 `frontend/.env` 中的 `VITE_API_URL` 配置
2. 检查后端是否启动：`curl http://localhost:8000/health`
3. 检查 CORS 配置：`backend/.env` 中的 `BACKEND_CORS_ORIGINS`

### MinIO 无法访问

```bash
# 检查MinIO服务
docker-compose ps minio

# 重启MinIO
docker-compose restart minio

# 检查MinIO日志
docker-compose logs -f minio
```

### Celery Worker 无法启动

```bash
# 检查Celery日志
docker-compose logs -f celery-worker

# 检查Redis连接
docker-compose exec redis redis-cli ping

# 重启Celery Worker
docker-compose restart celery-worker
```

---

## 🔒 安全配置

### 生产环境检查清单

在生产环境部署前，请确保：

- [ ] 修改所有默认密码
- [ ] 更新 `SECRET_KEY` 和 `JWT_SECRET_KEY`
- [ ] 设置 `DEBUG=False`
- [ ] 配置 HTTPS/SSL 证书
- [ ] 限制 CORS 允许的域名
- [ ] 配置防火墙规则
- [ ] 设置数据库备份策略
- [ ] 配置日志收集和监控
- [ ] 启用 Rate Limiting
- [ ] 配置文件上传大小限制

### 修改默认凭据

```bash
# PostgreSQL密码
# 修改 docker-compose.yml 中的 POSTGRES_PASSWORD

# MinIO凭据
# 修改 docker-compose.yml 中的 MINIO_ROOT_USER 和 MINIO_ROOT_PASSWORD

# JWT密钥
# 修改 backend/.env 中的 SECRET_KEY 和 JWT_SECRET_KEY

# 重新部署
docker-compose down
docker-compose up -d
```

---

## 📊 性能优化

### 后端性能

```bash
# 增加Celery Worker数量
docker-compose up -d --scale celery-worker=4

# 调整数据库连接池
# 修改 backend/.env 中的 DATABASE_POOL_SIZE
```

### 前端性能

```bash
# 生产构建
cd frontend
npm run build

# 使用Nginx部署静态文件
docker run -d -p 80:80 \
  -v $(pwd)/dist:/usr/share/nginx/html \
  nginx:alpine
```

---

## 📚 相关文档

- [API 文档](http://localhost:8000/api/v1/docs)
- [项目进度报告](PROJECT_PROGRESS_SUMMARY.md)
- [Phase 1 完成报告](PHASE1_COMPLETION_REPORT.md)
- [Phase 2 完成报告](PHASE2_COMPLETION_REPORT.md)
- [Phase 3 完成报告](PHASE3_COMPLETION_REPORT.md)

---

## 🤝 支持

遇到问题？

1. 查看[故障排除](#故障排除)章节
2. 查看日志：`docker-compose logs -f`
3. 提交 Issue
4. 查阅 API 文档

---

**更新时间**: 2025-11-21
**版本**: 1.0.0
