# NGSmodule 部署指南

本指南提供 NGSmodule 生产环境部署的完整步骤和最佳实践。

## 📋 目录

- [系统要求](#系统要求)
- [快速开始](#快速开始)
- [生产环境部署](#生产环境部署)
- [环境变量配置](#环境变量配置)
- [数据库设置](#数据库设置)
- [SSL/HTTPS 配置](#sslhttps-配置)
- [扩展和优化](#扩展和优化)
- [监控和日志](#监控和日志)
- [备份和恢复](#备份和恢复)
- [故障排除](#故障排除)

---

## 系统要求

### 最小要求
- **CPU**: 2 核心
- **内存**: 4 GB RAM
- **存储**: 20 GB 可用空间
- **操作系统**: Linux (Ubuntu 20.04+, CentOS 8+, Debian 11+)
- **Docker**: 20.10+
- **Docker Compose**: 2.0+

### 推荐配置
- **CPU**: 4+ 核心
- **内存**: 8+ GB RAM
- **存储**: 100+ GB SSD
- **网络**: 100 Mbps+

---

## 快速开始

### 1. 克隆仓库

```bash
git clone https://github.com/yourusername/NGSmodule.git
cd NGSmodule
```

### 2. 配置环境变量

```bash
# 复制环境变量模板
cp .env.example .env

# 编辑配置文件
nano .env
```

**必须修改的关键配置:**
```env
# 强密码 (使用: openssl rand -hex 32)
POSTGRES_PASSWORD=your_secure_postgres_password
REDIS_PASSWORD=your_secure_redis_password
SECRET_KEY=your_secret_key_32_chars_minimum
JWT_SECRET_KEY=your_jwt_secret_key_32_chars_minimum

# 生产环境域名
BACKEND_CORS_ORIGINS=https://yourdomain.com
VITE_API_URL=https://api.yourdomain.com/api/v1
VITE_WS_URL=wss://api.yourdomain.com/api/v1/ws
```

### 3. 启动服务

```bash
# 开发环境
docker-compose up -d

# 生产环境 (包含 Nginx)
docker-compose --profile production up -d
```

### 4. 初始化数据库

```bash
# 运行数据库迁移
docker-compose exec backend alembic upgrade head

# 创建初始管理员用户 (可选)
docker-compose exec backend python scripts/create_admin.py
```

### 5. 验证部署

```bash
# 检查所有服务状态
docker-compose ps

# 查看服务日志
docker-compose logs -f backend
```

访问: `http://localhost` (生产) 或 `http://localhost:3000` (开发)

---

## 生产环境部署

### 架构概览

```
                    ┌──────────────┐
                    │   Internet   │
                    └──────┬───────┘
                           │
                    ┌──────▼───────┐
                    │ Nginx (443)  │
                    │  SSL/TLS     │
                    └──┬────────┬──┘
                       │        │
           ┌───────────▼──┐  ┌──▼───────────┐
           │  Frontend    │  │  Backend     │
           │  (Nginx)     │  │  (FastAPI)   │
           └──────────────┘  └──┬────────┬──┘
                                │        │
                    ┌───────────▼──┐  ┌──▼──────────┐
                    │  PostgreSQL  │  │   Redis     │
                    └──────────────┘  └─────────────┘
```

### 步骤 1: 服务器准备

```bash
# 更新系统
sudo apt update && sudo apt upgrade -y

# 安装 Docker
curl -fsSL https://get.docker.com | sh
sudo usermod -aG docker $USER

# 安装 Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# 重新登录以应用组变更
exit
```

### 步骤 2: 防火墙配置

```bash
# 允许 HTTP/HTTPS
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# 允许 SSH (确保先允许 SSH!)
sudo ufw allow 22/tcp

# 启用防火墙
sudo ufw enable
```

### 步骤 3: 生产配置

编辑 `.env` 文件:

```env
# 生产模式
ENVIRONMENT=production
DEBUG=false
LOG_LEVEL=INFO

# 性能优化
UVICORN_WORKERS=4
CELERY_WORKERS=2
CELERY_CONCURRENCY=4

# 安全设置
BACKEND_CORS_ORIGINS=https://yourdomain.com
```

### 步骤 4: 启动生产环境

```bash
# 构建并启动所有服务
docker-compose --profile production build
docker-compose --profile production up -d

# 验证服务健康
docker-compose exec backend curl http://localhost:8000/api/v1/health
```

---

## 环境变量配置

### 核心配置说明

| 变量 | 说明 | 默认值 | 必填 |
|------|------|--------|------|
| `POSTGRES_PASSWORD` | PostgreSQL 密码 | - | ✅ |
| `REDIS_PASSWORD` | Redis 密码 | - | ✅ |
| `SECRET_KEY` | 应用密钥 | - | ✅ |
| `JWT_SECRET_KEY` | JWT 签名密钥 | - | ✅ |
| `BACKEND_CORS_ORIGINS` | 允许的前端域名 | - | ✅ |
| `VITE_API_URL` | 前端 API 地址 | - | ✅ |

### 生成安全密钥

```bash
# 生成 SECRET_KEY
openssl rand -hex 32

# 生成 JWT_SECRET_KEY
openssl rand -hex 32

# 生成强密码
openssl rand -base64 32
```

---

## 数据库设置

### 初始化数据库

数据库会自动初始化，包含:
- 扩展: `uuid-ossp`, `pg_trgm`
- 自定义类型: `project_status`, `task_status`
- 时区设置为 UTC

### 运行迁移

```bash
# 升级到最新版本
docker-compose exec backend alembic upgrade head

# 查看迁移历史
docker-compose exec backend alembic history

# 回滚迁移
docker-compose exec backend alembic downgrade -1
```

### 创建管理员用户

创建 `backend/scripts/create_admin.py`:

```python
import sys
sys.path.append('/app')

from app.db.session import SessionLocal
from app.models.user import User
from app.core.security import get_password_hash

db = SessionLocal()
admin = User(
    username="admin",
    email="admin@example.com",
    hashed_password=get_password_hash("admin123"),
    is_active=True,
    is_superuser=True
)
db.add(admin)
db.commit()
print("Admin user created successfully!")
```

运行:
```bash
docker-compose exec backend python scripts/create_admin.py
```

---

## SSL/HTTPS 配置

### 方案 1: Let's Encrypt (推荐)

使用 Certbot 自动获取免费 SSL 证书:

```bash
# 安装 Certbot
sudo apt install certbot python3-certbot-nginx -y

# 停止 Nginx 容器
docker-compose stop nginx

# 获取证书
sudo certbot certonly --standalone -d yourdomain.com -d www.yourdomain.com

# 证书位置
# /etc/letsencrypt/live/yourdomain.com/fullchain.pem
# /etc/letsencrypt/live/yourdomain.com/privkey.pem
```

更新 `docker-compose.yml` 中的 Nginx 配置:

```yaml
nginx:
  volumes:
    - /etc/letsencrypt:/etc/letsencrypt:ro
    - ./nginx/nginx-ssl.conf:/etc/nginx/nginx.conf:ro
```

创建 `nginx/nginx-ssl.conf`:

```nginx
server {
    listen 80;
    server_name yourdomain.com;
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl http2;
    server_name yourdomain.com;

    ssl_certificate /etc/letsencrypt/live/yourdomain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/yourdomain.com/privkey.pem;
    
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers HIGH:!aNULL:!MD5;
    ssl_prefer_server_ciphers on;

    # 前端
    location / {
        proxy_pass http://frontend:80;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    # 后端 API
    location /api/ {
        proxy_pass http://backend:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # WebSocket
    location /api/v1/ws {
        proxy_pass http://backend:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

### 方案 2: 自签名证书 (开发/测试)

```bash
# 生成自签名证书
mkdir -p nginx/ssl
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
  -keyout nginx/ssl/nginx-selfsigned.key \
  -out nginx/ssl/nginx-selfsigned.crt \
  -subj "/C=US/ST=State/L=City/O=Organization/CN=localhost"

# 挂载到容器
# docker-compose.yml:
# volumes:
#   - ./nginx/ssl:/etc/nginx/ssl:ro
```

### 自动续期

```bash
# 添加 crontab 任务
sudo crontab -e

# 每月 1 号凌晨 2 点续期
0 2 1 * * certbot renew --quiet && docker-compose restart nginx
```

---

## 扩展和优化

### 水平扩展

#### 扩展 Backend Workers

```bash
# 临时扩展到 4 个实例
docker-compose up -d --scale backend=4

# 配置负载均衡 (Nginx)
upstream backend_cluster {
    server backend_1:8000;
    server backend_2:8000;
    server backend_3:8000;
    server backend_4:8000;
}
```

#### 扩展 Celery Workers

```yaml
# docker-compose.yml
celery-worker:
  deploy:
    replicas: 4
```

### 性能优化

#### 数据库连接池

```env
# .env
DB_POOL_SIZE=20
DB_MAX_OVERFLOW=10
```

#### Redis 缓存配置

```python
# backend/app/core/config.py
REDIS_CACHE_TTL = 300  # 5 分钟
REDIS_MAX_CONNECTIONS = 50
```

#### Nginx 缓存

```nginx
# nginx.conf
proxy_cache_path /var/cache/nginx levels=1:2 keys_zone=api_cache:10m max_size=1g inactive=60m;

location /api/v1/public/ {
    proxy_cache api_cache;
    proxy_cache_valid 200 5m;
    proxy_pass http://backend:8000;
}
```

### 资源限制

```yaml
# docker-compose.yml
services:
  backend:
    deploy:
      resources:
        limits:
          cpus: '2'
          memory: 2G
        reservations:
          cpus: '1'
          memory: 1G
```

---

## 监控和日志

### 健康检查

所有服务都配置了健康检查:

```bash
# 检查服务健康状态
docker-compose ps

# Backend 健康端点
curl http://localhost:8000/api/v1/health

# PostgreSQL
docker-compose exec postgres pg_isready

# Redis
docker-compose exec redis redis-cli ping
```

### 日志管理

```bash
# 查看所有日志
docker-compose logs -f

# 查看特定服务
docker-compose logs -f backend

# 最近 100 行
docker-compose logs --tail=100 backend

# 导出日志
docker-compose logs --no-color > logs_$(date +%Y%m%d).txt
```

### 日志轮转

创建 `/etc/logrotate.d/docker-ngsmodule`:

```
/var/lib/docker/containers/*/*.log {
    rotate 7
    daily
    compress
    missingok
    delaycompress
    copytruncate
}
```

---

## 备份和恢复

### 自动备份

备份已通过 `celery-beat` 自动配置 (每日凌晨 2 点):

```bash
# 查看备份任务状态
docker-compose logs celery-beat | grep backup

# 手动触发备份
docker-compose exec backend bash /app/scripts/backup.sh
```

### 备份位置

```bash
# 默认备份目录
ls -lh /backups/

# 备份命名格式
# ngsmodule_backup_20250122_020000.sql.gz
```

### 恢复数据库

```bash
# 列出可用备份
docker-compose exec backend bash /app/scripts/restore.sh

# 恢复特定备份
docker-compose exec backend bash /app/scripts/restore.sh /backups/ngsmodule_backup_20250122_020000.sql.gz
```

### 手动备份

```bash
# 完整备份
docker-compose exec postgres pg_dump -U ngsmodule ngsmodule | gzip > backup_$(date +%Y%m%d).sql.gz

# 仅数据
docker-compose exec postgres pg_dump -U ngsmodule --data-only ngsmodule > data_backup.sql

# 仅结构
docker-compose exec postgres pg_dump -U ngsmodule --schema-only ngsmodule > schema_backup.sql
```

### 备份存储卷

```bash
# 备份 PostgreSQL 数据
docker run --rm \
  -v ngsmodule_postgres_data:/data \
  -v $(pwd):/backup \
  alpine tar czf /backup/postgres_volume_$(date +%Y%m%d).tar.gz -C /data .

# 备份应用存储
docker run --rm \
  -v ngsmodule_storage:/data \
  -v $(pwd):/backup \
  alpine tar czf /backup/storage_volume_$(date +%Y%m%d).tar.gz -C /data .
```

---

## 故障排除

### 常见问题

#### 1. 容器无法启动

```bash
# 检查日志
docker-compose logs backend

# 常见原因:
# - 端口被占用: sudo lsof -i :8000
# - 环境变量缺失: 检查 .env 文件
# - 权限问题: sudo chown -R $USER:$USER .
```

#### 2. 数据库连接失败

```bash
# 检查 PostgreSQL 状态
docker-compose exec postgres pg_isready

# 测试连接
docker-compose exec backend python -c "
from app.db.session import SessionLocal
db = SessionLocal()
print('Database connected!')
db.close()
"

# 检查密码
docker-compose exec postgres psql -U ngsmodule -d ngsmodule -c "SELECT 1"
```

#### 3. Frontend 无法访问 API

```bash
# 检查 CORS 配置
docker-compose logs backend | grep CORS

# 验证环境变量
docker-compose exec backend env | grep BACKEND_CORS_ORIGINS

# 测试 API 直接访问
curl http://localhost:8000/api/v1/health
```

#### 4. Celery 任务不执行

```bash
# 检查 Celery worker
docker-compose logs celery-worker

# 检查 Redis 连接
docker-compose exec backend python -c "
from app.core.celery_app import celery_app
print(celery_app.control.inspect().active())
"

# 手动执行任务测试
docker-compose exec backend celery -A app.core.celery_app worker --loglevel=debug
```

#### 5. 磁盘空间不足

```bash
# 清理 Docker 缓存
docker system prune -a --volumes

# 清理旧日志
docker-compose logs --no-color > /dev/null
truncate -s 0 /var/lib/docker/containers/*/*.log

# 清理旧备份
find /backups -name "ngsmodule_backup_*.sql.gz" -mtime +30 -delete
```

### 性能问题诊断

#### 慢查询分析

```sql
-- 启用慢查询日志 (PostgreSQL)
ALTER SYSTEM SET log_min_duration_statement = 1000; -- 记录超过 1 秒的查询
SELECT pg_reload_conf();

-- 查看慢查询
SELECT * FROM pg_stat_statements ORDER BY mean_exec_time DESC LIMIT 10;
```

#### 资源使用监控

```bash
# 容器资源使用
docker stats

# 数据库连接数
docker-compose exec postgres psql -U ngsmodule -d ngsmodule -c "
SELECT count(*) FROM pg_stat_activity;
"

# Redis 内存使用
docker-compose exec redis redis-cli INFO memory
```

### 紧急恢复流程

```bash
# 1. 停止所有服务
docker-compose down

# 2. 恢复最新备份
docker-compose up -d postgres
docker-compose exec backend bash /app/scripts/restore.sh /backups/latest_backup.sql.gz

# 3. 重启所有服务
docker-compose up -d

# 4. 验证数据完整性
docker-compose exec backend python scripts/verify_data.py
```

---

## 安全检查清单

- [ ] 修改所有默认密码
- [ ] 启用 HTTPS/SSL
- [ ] 配置防火墙规则
- [ ] 限制数据库外部访问
- [ ] 启用自动备份
- [ ] 配置日志轮转
- [ ] 更新 CORS 白名单
- [ ] 禁用 DEBUG 模式
- [ ] 设置资源限制
- [ ] 配置监控告警

---

## 维护计划

### 每日
- 检查备份状态
- 查看错误日志
- 监控磁盘空间

### 每周
- 审查访问日志
- 检查性能指标
- 更新依赖包

### 每月
- 测试备份恢复
- 安全补丁更新
- 容量规划评估

### 每季度
- 全面安全审计
- 灾难恢复演练
- 性能优化评估

---

## 附录

### 有用的命令

```bash
# 完全重置环境
docker-compose down -v
docker system prune -a

# 重建特定服务
docker-compose up -d --build backend

# 进入容器 shell
docker-compose exec backend bash
docker-compose exec postgres psql -U ngsmodule

# 导出/导入数据
docker-compose exec postgres pg_dump -U ngsmodule ngsmodule > dump.sql
docker-compose exec -T postgres psql -U ngsmodule ngsmodule < dump.sql
```

### 联系支持

- **问题反馈**: https://github.com/yourusername/NGSmodule/issues
- **文档**: https://docs.ngsmodule.com
- **邮件**: support@ngsmodule.com

---

**最后更新**: 2025-01-22  
**版本**: 1.0.0
