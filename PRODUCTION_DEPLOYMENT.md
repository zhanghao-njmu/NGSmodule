# NGSmodule 生产环境部署指南

## 📋 目录

- [系统要求](#系统要求)
- [快速部署](#快速部署)
- [详细部署步骤](#详细部署步骤)
- [配置说明](#配置说明)
- [运维管理](#运维管理)
- [监控和日志](#监控和日志)
- [备份和恢复](#备份和恢复)
- [故障排查](#故障排查)
- [性能优化](#性能优化)
- [安全加固](#安全加固)

---

## 系统要求

### 硬件要求

**最低配置：**
- CPU: 4核心
- 内存: 8GB RAM
- 存储: 100GB SSD
- 网络: 100Mbps

**推荐配置：**
- CPU: 8核心或更多
- 内存: 16GB RAM或更多
- 存储: 500GB SSD或更多
- 网络: 1Gbps

### 软件要求

- **操作系统**: Ubuntu 20.04+ / CentOS 8+ / Debian 11+
- **Docker**: 20.10.0 或更高
- **Docker Compose**: 2.0.0 或更高
- **Git**: 2.25.0 或更高
- **OpenSSL**: 1.1.1 或更高

### 网络要求

- 公网IP地址
- 域名（推荐）
- 开放端口: 80 (HTTP), 443 (HTTPS)
- 防火墙配置正确

---

## 快速部署

### 一键部署（推荐新用户）

```bash
# 1. 克隆代码
git clone https://github.com/yourusername/NGSmodule.git
cd NGSmodule

# 2. 运行自动化部署脚本
./deploy-production.sh setup
./deploy-production.sh build
./deploy-production.sh start

# 3. 检查状态
./deploy-production.sh status
```

### 手动部署

```bash
# 1. 创建环境配置
cp .env.production .env

# 2. 编辑配置文件（必需！）
nano .env  # 更新所有 CHANGE_THIS 占位符

# 3. 创建必需目录
mkdir -p backups/postgres ssl logs
chmod 700 backups ssl logs

# 4. 生成密钥
openssl rand -hex 32  # 用于 SECRET_KEY
openssl rand -hex 32  # 用于 JWT_SECRET_KEY

# 5. 构建和启动
docker-compose -f docker-compose.prod.yml build
docker-compose -f docker-compose.prod.yml up -d

# 6. 验证部署
curl http://localhost/api/health
```

---

## 详细部署步骤

### 步骤 1: 准备服务器

```bash
# 更新系统
sudo apt update && sudo apt upgrade -y

# 安装必需软件
sudo apt install -y git curl wget openssl

# 安装 Docker
curl -fsSL https://get.docker.com | sh
sudo usermod -aG docker $USER

# 安装 Docker Compose
sudo curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose

# 验证安装
docker --version
docker-compose --version
```

### 步骤 2: 配置防火墙

```bash
# UFW (Ubuntu/Debian)
sudo ufw allow 22/tcp    # SSH
sudo ufw allow 80/tcp    # HTTP
sudo ufw allow 443/tcp   # HTTPS
sudo ufw enable

# 或 firewalld (CentOS/RHEL)
sudo firewall-cmd --permanent --add-service=http
sudo firewall-cmd --permanent --add-service=https
sudo firewall-cmd --reload
```

### 步骤 3: 配置 SSL 证书

#### 选项 A: Let's Encrypt（推荐）

```bash
# 安装 Certbot
sudo apt install -y certbot

# 获取证书（需要域名指向服务器）
sudo certbot certonly --standalone -d yourdomain.com -d www.yourdomain.com

# 复制证书到项目目录
sudo cp /etc/letsencrypt/live/yourdomain.com/fullchain.pem ./ssl/cert.pem
sudo cp /etc/letsencrypt/live/yourdomain.com/privkey.pem ./ssl/key.pem
sudo chown $USER:$USER ./ssl/*.pem
```

#### 选项 B: 自签名证书（测试用）

```bash
mkdir -p ssl
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
  -keyout ssl/key.pem -out ssl/cert.pem \
  -subj "/C=US/ST=State/L=City/O=Organization/CN=yourdomain.com"
```

### 步骤 4: 配置环境变量

```bash
# 复制模板
cp .env.production .env

# 生成安全密钥
echo "SECRET_KEY=$(openssl rand -hex 32)" >> .env.generated
echo "JWT_SECRET_KEY=$(openssl rand -hex 32)" >> .env.generated
echo "POSTGRES_PASSWORD=$(openssl rand -base64 32)" >> .env.generated
echo "REDIS_PASSWORD=$(openssl rand -base64 32)" >> .env.generated

# 编辑 .env 文件
nano .env
```

**必须修改的关键配置：**

```bash
# 域名配置
VITE_API_URL=https://yourdomain.com/api/v1
VITE_WS_URL=wss://yourdomain.com/ws

# CORS配置
BACKEND_CORS_ORIGINS=https://yourdomain.com,https://www.yourdomain.com

# 数据库密码
POSTGRES_PASSWORD=<生成的强密码>

# Redis密码
REDIS_PASSWORD=<生成的强密码>

# MinIO配置
MINIO_ROOT_USER=admin
MINIO_ROOT_PASSWORD=<生成的强密码>

# 应用密钥
SECRET_KEY=<生成的64位hex>
JWT_SECRET_KEY=<生成的64位hex>
```

### 步骤 5: 构建和启动

```bash
# 构建 Docker 镜像
docker-compose -f docker-compose.prod.yml build --no-cache

# 启动所有服务
docker-compose -f docker-compose.prod.yml up -d

# 查看日志
docker-compose -f docker-compose.prod.yml logs -f

# 检查服务状态
docker-compose -f docker-compose.prod.yml ps
```

### 步骤 6: 初始化数据库

```bash
# 进入后端容器
docker-compose -f docker-compose.prod.yml exec backend bash

# 运行数据库迁移（如果有）
# python manage.py migrate

# 创建管理员用户（如果需要）
# python manage.py createsuperuser

# 退出容器
exit
```

### 步骤 7: 验证部署

```bash
# 检查所有服务健康状态
./deploy-production.sh status

# 测试 API
curl http://localhost/api/health
curl https://yourdomain.com/api/health

# 测试前端
curl http://localhost/
curl https://yourdomain.com/

# 检查数据库连接
docker-compose -f docker-compose.prod.yml exec postgres psql -U ngsmodule -d ngsmodule -c "SELECT version();"

# 检查 Redis
docker-compose -f docker-compose.prod.yml exec redis redis-cli -a $REDIS_PASSWORD ping
```

---

## 配置说明

### Nginx 反向代理

主要配置文件：`nginx-proxy.conf`

**关键配置：**

```nginx
# 文件上传大小限制
client_max_body_size 100M;

# WebSocket支持
location /ws/ {
    proxy_pass http://backend_api;
    proxy_http_version 1.1;
    proxy_set_header Upgrade $http_upgrade;
    proxy_set_header Connection "upgrade";
}

# API速率限制
limit_req zone=api_limit burst=20 nodelay;

# 认证端点严格限制
limit_req zone=auth_limit burst=5 nodelay;
```

### Docker Compose 服务

生产配置文件：`docker-compose.prod.yml`

**服务架构：**

```
nginx (反向代理)
  ├── frontend (React SPA)
  └── backend (FastAPI)
       ├── postgres (数据库)
       ├── redis (缓存/队列)
       ├── minio (对象存储)
       ├── celery-worker (异步任务)
       ├── celery-beat (定时任务)
       └── flower (监控)
```

### 环境变量

完整配置见 `.env.production`

**关键变量分类：**

1. **安全配置**
   - SECRET_KEY
   - JWT_SECRET_KEY
   - 各种 PASSWORD

2. **数据库配置**
   - POSTGRES_*
   - REDIS_*
   - MINIO_*

3. **应用配置**
   - ENVIRONMENT=production
   - DEBUG=false
   - LOG_LEVEL=INFO

4. **域名配置**
   - VITE_API_URL
   - VITE_WS_URL
   - BACKEND_CORS_ORIGINS

---

## 运维管理

### 日常运维命令

```bash
# 查看所有服务状态
docker-compose -f docker-compose.prod.yml ps

# 查看实时日志
docker-compose -f docker-compose.prod.yml logs -f [service_name]

# 重启特定服务
docker-compose -f docker-compose.prod.yml restart backend

# 重启所有服务
docker-compose -f docker-compose.prod.yml restart

# 停止所有服务
docker-compose -f docker-compose.prod.yml down

# 停止并删除数据卷（危险！）
docker-compose -f docker-compose.prod.yml down -v
```

### 服务扩容

```bash
# 水平扩展 Celery workers
docker-compose -f docker-compose.prod.yml up -d --scale celery-worker=4

# 增加后端实例（需要负载均衡配置）
docker-compose -f docker-compose.prod.yml up -d --scale backend=3
```

### 更新部署

```bash
# 使用自动化脚本
./deploy-production.sh update

# 或手动更新
git pull
docker-compose -f docker-compose.prod.yml build
docker-compose -f docker-compose.prod.yml up -d --no-deps backend frontend
```

---

## 监控和日志

### 日志管理

**查看日志：**

```bash
# 所有服务日志
docker-compose -f docker-compose.prod.yml logs -f

# 特定服务日志
docker-compose -f docker-compose.prod.yml logs -f backend

# 最近100行日志
docker-compose -f docker-compose.prod.yml logs --tail=100 backend

# 导出日志
docker-compose -f docker-compose.prod.yml logs --no-color > logs/docker-$(date +%Y%m%d).log
```

**日志配置：**

所有服务使用 json-file 驱动，配置了日志轮转：
- max-size: 10-50MB
- max-file: 3-5个文件

### Celery 监控（Flower）

访问：http://localhost:5555

默认凭据在 `.env` 中配置：
- 用户名: FLOWER_USER
- 密码: FLOWER_PASSWORD

### Nginx 状态监控

```bash
# 访问 Nginx 状态页
curl http://localhost:8080/nginx-status
```

### 健康检查

```bash
# API健康检查
curl http://localhost/api/health

# 数据库健康检查
docker-compose -f docker-compose.prod.yml exec postgres pg_isready

# Redis健康检查
docker-compose -f docker-compose.prod.yml exec redis redis-cli ping

# 所有服务健康状态
./deploy-production.sh status
```

---

## 备份和恢复

### 自动备份

**配置 cron 任务：**

```bash
# 编辑 crontab
crontab -e

# 添加每日凌晨2点备份
0 2 * * * cd /path/to/NGSmodule && ./deploy-production.sh backup >> logs/backup.log 2>&1
```

### 手动备份

```bash
# 使用脚本备份
./deploy-production.sh backup

# 或手动备份
docker-compose -f docker-compose.prod.yml exec -T postgres \
  pg_dump -U ngsmodule ngsmodule | gzip > backups/postgres/backup_$(date +%Y%m%d_%H%M%S).sql.gz
```

### 数据恢复

```bash
# 使用脚本恢复
./deploy-production.sh restore backups/postgres/backup_20240101_020000.sql.gz

# 或手动恢复
gunzip -c backups/postgres/backup_20240101_020000.sql.gz | \
  docker-compose -f docker-compose.prod.yml exec -T postgres \
  psql -U ngsmodule ngsmodule
```

### 完整备份策略

**需要备份的内容：**

1. **数据库**（PostgreSQL）
   - 位置: backups/postgres/
   - 频率: 每日
   - 保留: 30天

2. **持久化数据卷**
   - postgres_data
   - redis_data
   - minio_data
   - upload_data
   - result_data

3. **配置文件**
   - .env（包含密钥，需加密保存）
   - nginx-proxy.conf
   - docker-compose.prod.yml

4. **SSL证书**
   - ssl/cert.pem
   - ssl/key.pem

**备份命令：**

```bash
# 备份所有数据卷
docker run --rm \
  -v ngsmodule_postgres_data:/data/postgres \
  -v ngsmodule_redis_data:/data/redis \
  -v ngsmodule_minio_data:/data/minio \
  -v $(pwd)/backups:/backup \
  alpine tar czf /backup/volumes_$(date +%Y%m%d).tar.gz /data
```

---

## 故障排查

### 常见问题

#### 1. 容器无法启动

```bash
# 查看容器日志
docker-compose -f docker-compose.prod.yml logs [service_name]

# 检查容器状态
docker-compose -f docker-compose.prod.yml ps

# 检查端口占用
sudo netstat -tulpn | grep :80
sudo netstat -tulpn | grep :8000
```

#### 2. 数据库连接失败

```bash
# 检查数据库健康状态
docker-compose -f docker-compose.prod.yml exec postgres pg_isready

# 测试数据库连接
docker-compose -f docker-compose.prod.yml exec postgres \
  psql -U ngsmodule -d ngsmodule -c "SELECT 1;"

# 检查环境变量
docker-compose -f docker-compose.prod.yml exec backend env | grep DATABASE
```

#### 3. 前端无法访问API

```bash
# 检查CORS配置
grep BACKEND_CORS_ORIGINS .env

# 检查前端环境变量
docker-compose -f docker-compose.prod.yml exec frontend env | grep VITE

# 检查Nginx配置
docker-compose -f docker-compose.prod.yml exec nginx nginx -t

# 重新加载Nginx
docker-compose -f docker-compose.prod.yml exec nginx nginx -s reload
```

#### 4. 上传文件失败

```bash
# 检查文件大小限制
grep client_max_body_size nginx-proxy.conf

# 检查存储卷
docker volume ls
docker volume inspect ngsmodule_upload_data

# 检查权限
docker-compose -f docker-compose.prod.yml exec backend ls -la /data/uploads
```

#### 5. SSL证书问题

```bash
# 检查证书有效期
openssl x509 -in ssl/cert.pem -noout -dates

# 测试SSL配置
openssl s_client -connect yourdomain.com:443 -servername yourdomain.com

# 更新证书（Let's Encrypt）
sudo certbot renew
```

### 诊断工具

```bash
# 系统资源使用
docker stats

# 容器详细信息
docker inspect [container_name]

# 网络诊断
docker network inspect ngsmodule_ngsmodule_network

# 进入容器调试
docker-compose -f docker-compose.prod.yml exec backend bash
```

---

## 性能优化

### 数据库优化

**PostgreSQL配置优化：**

```sql
-- 查看当前配置
SHOW all;

-- 优化查询
EXPLAIN ANALYZE SELECT * FROM ...;

-- 创建索引
CREATE INDEX idx_name ON table_name(column_name);

-- 清理数据库
VACUUM ANALYZE;
```

### Redis优化

```bash
# 内存使用情况
docker-compose -f docker-compose.prod.yml exec redis redis-cli info memory

# 清理过期键
docker-compose -f docker-compose.prod.yml exec redis redis-cli --scan --pattern "*" | xargs redis-cli del

# 优化配置（在docker-compose.prod.yml中）
command: >
  redis-server
  --maxmemory 512mb
  --maxmemory-policy allkeys-lru
```

### Nginx优化

**worker进程调优：**

```nginx
# 在nginx-proxy.conf添加
worker_processes auto;
worker_rlimit_nofile 65535;

events {
    worker_connections 4096;
    use epoll;
    multi_accept on;
}

http {
    # 启用HTTP/2
    listen 443 ssl http2;

    # 连接池
    keepalive_timeout 65;
    keepalive_requests 100;

    # 缓冲区
    client_body_buffer_size 128k;
    client_header_buffer_size 1k;
    large_client_header_buffers 4 16k;
}
```

### 应用层优化

**后端优化：**

```python
# 数据库连接池
DB_POOL_SIZE=20
DB_MAX_OVERFLOW=10

# Uvicorn workers
UVICORN_WORKERS=4  # CPU核心数

# Celery并发
CELERY_CONCURRENCY=4
```

**前端优化：**

已实现的优化：
- ✅ 路由懒加载（-55% bundle）
- ✅ Gzip压缩
- ✅ 静态资源缓存（1年）
- ✅ 代码分割（14个chunks）

---

## 安全加固

### 系统安全

```bash
# 配置自动安全更新
sudo apt install unattended-upgrades
sudo dpkg-reconfigure -plow unattended-upgrades

# 配置fail2ban
sudo apt install fail2ban
sudo systemctl enable fail2ban

# SSH安全加固
sudo nano /etc/ssh/sshd_config
# 设置：
# PermitRootLogin no
# PasswordAuthentication no
# PubkeyAuthentication yes
```

### 应用安全

**1. 定期更新依赖：**

```bash
# 更新前端依赖
cd frontend && npm audit fix

# 更新后端依赖
cd backend && pip list --outdated
```

**2. 容器安全扫描：**

```bash
# 使用docker scout（如果可用）
docker scout cve ngsmodule/backend:latest
docker scout cve ngsmodule/frontend:latest
```

**3. 密钥管理：**

```bash
# .env文件权限
chmod 600 .env

# 使用环境变量注入（推荐生产环境）
# 不要在代码中硬编码密钥
```

**4. 网络安全：**

```bash
# 限制数据库外部访问
# 在docker-compose.prod.yml中，移除postgres的ports映射

# 使用内部网络
# 只有nginx暴露80/443端口
```

### 安全清单

- [ ] 所有密码使用强随机生成
- [ ] .env文件权限设置为600
- [ ] SSL证书配置正确
- [ ] 防火墙只开放必要端口
- [ ] 数据库不对外暴露
- [ ] Redis需要密码认证
- [ ] 启用CORS限制
- [ ] 配置CSP头
- [ ] 启用HSTS
- [ ] 定期备份测试
- [ ] 日志监控告警
- [ ] 定期安全审计

---

## 维护计划

### 每日检查

```bash
# 自动化健康检查脚本
./deploy-production.sh status

# 检查磁盘空间
df -h

# 检查日志错误
docker-compose -f docker-compose.prod.yml logs --since 24h | grep -i error
```

### 每周维护

- 检查备份完整性
- 审查系统日志
- 更新系统补丁
- 检查性能指标

### 每月维护

- 清理旧日志文件
- 优化数据库
- 更新依赖包
- 安全审计

---

## 联系和支持

- **文档**: https://github.com/yourusername/NGSmodule/wiki
- **问题跟踪**: https://github.com/yourusername/NGSmodule/issues
- **邮件**: admin@yourdomain.com

---

**最后更新**: 2024-01-01
**版本**: 1.0.0
**维护者**: NGSmodule Team
