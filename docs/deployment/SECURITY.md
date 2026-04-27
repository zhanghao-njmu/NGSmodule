# NGSmodule 安全最佳实践指南

本文档提供 NGSmodule 生产环境的安全加固建议和最佳实践。

## 📋 目录

- [认证和授权](#认证和授权)
- [密码策略](#密码策略)
- [API 安全](#api-安全)
- [数据库安全](#数据库安全)
- [文件上传安全](#文件上传安全)
- [网络安全](#网络安全)
- [容器安全](#容器安全)
- [密钥管理](#密钥管理)
- [安全监控](#安全监控)
- [安全检查清单](#安全检查清单)

---

## 认证和授权

### JWT Token 安全

#### 1. Token 配置

```env
# .env - 强密钥配置
JWT_SECRET_KEY=使用至少32字符的随机密钥
JWT_ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30  # 短期访问令牌
REFRESH_TOKEN_EXPIRE_DAYS=7     # 刷新令牌
```

生成安全密钥:
```bash
openssl rand -hex 32
```

#### 2. Token 刷新机制

```python
# backend/app/api/v1/auth.py
@router.post("/refresh")
async def refresh_token(
    refresh_token: str = Body(...),
    db: Session = Depends(get_db)
):
    # 验证刷新令牌
    payload = verify_token(refresh_token)
    if payload.get("type") != "refresh":
        raise HTTPException(status_code=401, detail="Invalid token type")
    
    # 检查令牌是否在黑名单中
    if await is_token_blacklisted(refresh_token):
        raise HTTPException(status_code=401, detail="Token revoked")
    
    # 生成新的访问令牌
    user = get_user(db, user_id=payload["sub"])
    access_token = create_access_token(user.id)
    
    return {"access_token": access_token, "token_type": "bearer"}
```

#### 3. Token 黑名单 (Redis)

```python
# backend/app/core/security.py
import redis
from datetime import timedelta

redis_client = redis.Redis(
    host=settings.REDIS_HOST,
    password=settings.REDIS_PASSWORD,
    decode_responses=True
)

async def blacklist_token(token: str, expires_in: int):
    """将 token 加入黑名单"""
    redis_client.setex(
        f"blacklist:{token}",
        timedelta(seconds=expires_in),
        "1"
    )

async def is_token_blacklisted(token: str) -> bool:
    """检查 token 是否在黑名单中"""
    return redis_client.exists(f"blacklist:{token}") > 0
```

### 权限控制

#### 1. 基于角色的访问控制 (RBAC)

```python
# backend/app/models/user.py
class User(Base):
    __tablename__ = "users"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    username = Column(String, unique=True, index=True, nullable=False)
    email = Column(String, unique=True, index=True, nullable=False)
    hashed_password = Column(String, nullable=False)
    is_active = Column(Boolean, default=True)
    is_superuser = Column(Boolean, default=False)
    role = Column(String, default="user")  # user, admin, viewer
```

#### 2. 权限装饰器

```python
# backend/app/core/deps.py
from functools import wraps
from fastapi import HTTPException, status

def require_role(*allowed_roles: str):
    """要求特定角色才能访问"""
    def decorator(func):
        @wraps(func)
        async def wrapper(*args, current_user: User = Depends(get_current_user), **kwargs):
            if current_user.role not in allowed_roles and not current_user.is_superuser:
                raise HTTPException(
                    status_code=status.HTTP_403_FORBIDDEN,
                    detail="Insufficient permissions"
                )
            return await func(*args, current_user=current_user, **kwargs)
        return wrapper
    return decorator

# 使用示例
@router.delete("/projects/{project_id}")
@require_role("admin")
async def delete_project(project_id: UUID, current_user: User = Depends(get_current_user)):
    pass
```

---

## 密码策略

### 1. 密码强度要求

```python
# backend/app/utils/validation.py
import re
from fastapi import HTTPException

def validate_password_strength(password: str):
    """
    验证密码强度:
    - 至少 8 个字符
    - 包含大写字母
    - 包含小写字母
    - 包含数字
    - 包含特殊字符
    """
    if len(password) < 8:
        raise HTTPException(
            status_code=400,
            detail="密码至少需要 8 个字符"
        )
    
    if not re.search(r"[A-Z]", password):
        raise HTTPException(
            status_code=400,
            detail="密码必须包含至少一个大写字母"
        )
    
    if not re.search(r"[a-z]", password):
        raise HTTPException(
            status_code=400,
            detail="密码必须包含至少一个小写字母"
        )
    
    if not re.search(r"\d", password):
        raise HTTPException(
            status_code=400,
            detail="密码必须包含至少一个数字"
        )
    
    if not re.search(r"[!@#$%^&*(),.?\":{}|<>]", password):
        raise HTTPException(
            status_code=400,
            detail="密码必须包含至少一个特殊字符"
        )
    
    return True
```

### 2. 密码哈希

```python
# backend/app/core/security.py
from passlib.context import CryptContext

# 使用 bcrypt 算法
pwd_context = CryptContext(
    schemes=["bcrypt"],
    deprecated="auto",
    bcrypt__rounds=12  # 增加计算成本
)

def get_password_hash(password: str) -> str:
    """生成密码哈希"""
    return pwd_context.hash(password)

def verify_password(plain_password: str, hashed_password: str) -> bool:
    """验证密码"""
    return pwd_context.verify(plain_password, hashed_password)
```

### 3. 登录尝试限制

```python
# backend/app/api/v1/auth.py
from datetime import datetime, timedelta

MAX_LOGIN_ATTEMPTS = 5
LOCKOUT_DURATION = timedelta(minutes=15)

async def check_login_attempts(username: str, redis_client) -> bool:
    """检查登录尝试次数"""
    key = f"login_attempts:{username}"
    attempts = redis_client.get(key)
    
    if attempts and int(attempts) >= MAX_LOGIN_ATTEMPTS:
        # 检查锁定时间
        ttl = redis_client.ttl(key)
        if ttl > 0:
            raise HTTPException(
                status_code=429,
                detail=f"Too many login attempts. Try again in {ttl} seconds."
            )
    
    return True

async def record_failed_login(username: str, redis_client):
    """记录失败的登录尝试"""
    key = f"login_attempts:{username}"
    redis_client.incr(key)
    redis_client.expire(key, int(LOCKOUT_DURATION.total_seconds()))

async def clear_login_attempts(username: str, redis_client):
    """清除登录尝试记录"""
    redis_client.delete(f"login_attempts:{username}")
```

---

## API 安全

### 1. 速率限制

```python
# backend/app/middleware/rate_limit.py
from fastapi import Request, HTTPException
from datetime import datetime, timedelta
import redis

class RateLimitMiddleware:
    def __init__(self, redis_client, max_requests: int = 100, window: int = 60):
        self.redis = redis_client
        self.max_requests = max_requests
        self.window = window
    
    async def __call__(self, request: Request, call_next):
        # 获取客户端 IP
        client_ip = request.client.host
        key = f"rate_limit:{client_ip}"
        
        # 检查请求次数
        current = self.redis.get(key)
        
        if current and int(current) >= self.max_requests:
            raise HTTPException(
                status_code=429,
                detail="Too many requests"
            )
        
        # 增加计数
        pipe = self.redis.pipeline()
        pipe.incr(key)
        pipe.expire(key, self.window)
        pipe.execute()
        
        response = await call_next(request)
        return response

# backend/app/main.py
from app.middleware.rate_limit import RateLimitMiddleware

app.add_middleware(
    RateLimitMiddleware,
    redis_client=redis_client,
    max_requests=100,
    window=60
)
```

### 2. CORS 配置

```python
# backend/app/main.py
from fastapi.middleware.cors import CORSMiddleware

# 严格的 CORS 配置
origins = [
    "https://yourdomain.com",
    "https://www.yourdomain.com",
]

# 开发环境可以添加 localhost
if settings.ENVIRONMENT == "development":
    origins.extend([
        "http://localhost:3000",
        "http://localhost:5173",
    ])

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE"],  # 明确指定方法
    allow_headers=["Content-Type", "Authorization"],  # 明确指定头部
    max_age=600,  # 预检请求缓存时间
)
```

### 3. 输入验证

```python
# backend/app/schemas/project.py
from pydantic import BaseModel, validator, constr
import re

class ProjectCreate(BaseModel):
    name: constr(min_length=1, max_length=100)  # 限制长度
    description: str | None = None
    
    @validator('name')
    def validate_name(cls, v):
        # 只允许字母、数字、下划线和连字符
        if not re.match(r'^[\w-]+$', v):
            raise ValueError('Project name can only contain letters, numbers, underscores and hyphens')
        return v
    
    @validator('description')
    def validate_description(cls, v):
        if v and len(v) > 1000:
            raise ValueError('Description too long (max 1000 characters)')
        return v

    class Config:
        # 严格模式：不允许额外字段
        extra = 'forbid'
```

### 4. SQL 注入防护

```python
# 使用 SQLAlchemy ORM (参数化查询)
# ✅ 安全
projects = db.query(Project).filter(Project.name == user_input).all()

# ❌ 不安全 - 永远不要这样做!
# db.execute(f"SELECT * FROM projects WHERE name = '{user_input}'")

# 如果必须使用原始 SQL，使用参数化查询
from sqlalchemy import text
stmt = text("SELECT * FROM projects WHERE name = :name")
results = db.execute(stmt, {"name": user_input})
```

---

## 数据库安全

### 1. 连接安全

```env
# .env
POSTGRES_USER=ngsmodule
POSTGRES_PASSWORD=使用强密码_至少16字符
POSTGRES_DB=ngsmodule
POSTGRES_HOST=postgres  # 仅容器内部访问

# 生产环境不要暴露数据库端口
# 移除 docker-compose.yml 中的:
# ports:
#   - "5432:5432"
```

### 2. 数据库用户权限

```sql
-- 创建只读用户 (用于报表等)
CREATE USER ngs_readonly WITH PASSWORD 'strong_password';
GRANT CONNECT ON DATABASE ngsmodule TO ngs_readonly;
GRANT USAGE ON SCHEMA public TO ngs_readonly;
GRANT SELECT ON ALL TABLES IN SCHEMA public TO ngs_readonly;

-- 创建应用用户 (限制权限)
CREATE USER ngs_app WITH PASSWORD 'strong_password';
GRANT CONNECT ON DATABASE ngsmodule TO ngs_app;
GRANT USAGE, CREATE ON SCHEMA public TO ngs_app;
GRANT SELECT, INSERT, UPDATE, DELETE ON ALL TABLES IN SCHEMA public TO ngs_app;

-- 撤销危险权限
REVOKE CREATE ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON DATABASE ngsmodule FROM PUBLIC;
```

### 3. 数据加密

```python
# backend/app/core/encryption.py
from cryptography.fernet import Fernet
import base64

class FieldEncryption:
    def __init__(self, key: str):
        self.cipher = Fernet(key.encode())
    
    def encrypt(self, data: str) -> str:
        """加密敏感字段"""
        return self.cipher.encrypt(data.encode()).decode()
    
    def decrypt(self, encrypted_data: str) -> str:
        """解密敏感字段"""
        return self.cipher.decrypt(encrypted_data.encode()).decode()

# 使用示例
encryptor = FieldEncryption(settings.ENCRYPTION_KEY)

# 加密敏感配置
user.api_key_encrypted = encryptor.encrypt(api_key)
```

### 4. 审计日志

```python
# backend/app/models/audit_log.py
class AuditLog(Base):
    __tablename__ = "audit_logs"
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.id"))
    action = Column(String, nullable=False)  # CREATE, UPDATE, DELETE
    resource_type = Column(String, nullable=False)  # Project, Sample, etc.
    resource_id = Column(UUID(as_uuid=True))
    changes = Column(JSON)  # 记录变更内容
    ip_address = Column(String)
    user_agent = Column(String)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

# 记录操作
async def log_audit(
    db: Session,
    user_id: UUID,
    action: str,
    resource_type: str,
    resource_id: UUID,
    changes: dict,
    request: Request
):
    audit = AuditLog(
        user_id=user_id,
        action=action,
        resource_type=resource_type,
        resource_id=resource_id,
        changes=changes,
        ip_address=request.client.host,
        user_agent=request.headers.get("user-agent")
    )
    db.add(audit)
    db.commit()
```

---

## 文件上传安全

### 1. 文件类型验证

```python
# backend/app/utils/file_validation.py
from fastapi import HTTPException, UploadFile
import magic
import os

ALLOWED_EXTENSIONS = {
    'fastq', 'fq', 'fastq.gz', 'fq.gz',
    'bam', 'sam', 'vcf', 'vcf.gz',
    'csv', 'txt'
}

ALLOWED_MIME_TYPES = {
    'application/gzip',
    'text/plain',
    'application/octet-stream',
}

async def validate_file_upload(file: UploadFile, max_size: int = 100 * 1024 * 1024):
    """
    验证上传文件
    - 检查文件扩展名
    - 检查 MIME 类型
    - 检查文件大小
    """
    # 检查文件扩展名
    filename = file.filename.lower()
    if not any(filename.endswith(ext) for ext in ALLOWED_EXTENSIONS):
        raise HTTPException(
            status_code=400,
            detail=f"File type not allowed. Allowed types: {', '.join(ALLOWED_EXTENSIONS)}"
        )
    
    # 读取文件头部检查 MIME 类型
    file_header = await file.read(2048)
    await file.seek(0)  # 重置文件指针
    
    mime_type = magic.from_buffer(file_header, mime=True)
    if mime_type not in ALLOWED_MIME_TYPES:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid file type. Detected: {mime_type}"
        )
    
    # 检查文件大小
    await file.seek(0, os.SEEK_END)
    file_size = await file.tell()
    await file.seek(0)
    
    if file_size > max_size:
        raise HTTPException(
            status_code=400,
            detail=f"File too large. Maximum size: {max_size / 1024 / 1024} MB"
        )
    
    return True
```

### 2. 安全的文件名处理

```python
# backend/app/utils/file_utils.py
import uuid
import os
import re

def sanitize_filename(filename: str) -> str:
    """清理文件名，防止路径遍历攻击"""
    # 移除路径组件
    filename = os.path.basename(filename)
    
    # 只保留字母、数字、点和下划线
    filename = re.sub(r'[^\w\.-]', '_', filename)
    
    # 限制长度
    name, ext = os.path.splitext(filename)
    if len(name) > 100:
        name = name[:100]
    
    return f"{name}{ext}"

def generate_safe_filepath(upload_dir: str, original_filename: str) -> str:
    """生成安全的文件路径"""
    # 清理文件名
    safe_filename = sanitize_filename(original_filename)
    
    # 使用 UUID 避免文件名冲突
    name, ext = os.path.splitext(safe_filename)
    unique_filename = f"{uuid.uuid4().hex}_{name}{ext}"
    
    # 确保上传目录存在且安全
    upload_dir = os.path.abspath(upload_dir)
    filepath = os.path.join(upload_dir, unique_filename)
    
    # 防止路径遍历
    if not filepath.startswith(upload_dir):
        raise ValueError("Invalid file path detected")
    
    return filepath
```

### 3. 病毒扫描 (可选)

```python
# backend/app/utils/virus_scan.py
import subprocess
from fastapi import HTTPException

async def scan_file_for_viruses(filepath: str):
    """使用 ClamAV 扫描文件"""
    try:
        result = subprocess.run(
            ['clamscan', '--no-summary', filepath],
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            # 发现病毒
            os.remove(filepath)
            raise HTTPException(
                status_code=400,
                detail="File rejected: potential security threat detected"
            )
    except subprocess.TimeoutExpired:
        raise HTTPException(
            status_code=500,
            detail="Virus scan timeout"
        )
    except FileNotFoundError:
        # ClamAV 未安装，记录警告但不阻止上传
        logger.warning("ClamAV not installed, skipping virus scan")
```

---

## 网络安全

### 1. HTTPS 强制

```nginx
# nginx/nginx.conf
server {
    listen 80;
    server_name yourdomain.com;
    
    # 重定向所有 HTTP 到 HTTPS
    return 301 https://$server_name$request_uri;
}

server {
    listen 443 ssl http2;
    server_name yourdomain.com;
    
    # SSL 配置
    ssl_certificate /etc/letsencrypt/live/yourdomain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/yourdomain.com/privkey.pem;
    
    # 强 SSL 配置
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers 'ECDHE-ECDSA-AES128-GCM-SHA256:ECDHE-RSA-AES128-GCM-SHA256:ECDHE-ECDSA-AES256-GCM-SHA384:ECDHE-RSA-AES256-GCM-SHA384';
    ssl_prefer_server_ciphers on;
    
    # HSTS
    add_header Strict-Transport-Security "max-age=31536000; includeSubDomains" always;
    
    # 其他安全头
    add_header X-Frame-Options "SAMEORIGIN" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;
    add_header Referrer-Policy "no-referrer-when-downgrade" always;
    add_header Content-Security-Policy "default-src 'self'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; style-src 'self' 'unsafe-inline';" always;
}
```

### 2. 防火墙规则

```bash
# UFW 防火墙配置
sudo ufw default deny incoming
sudo ufw default allow outgoing

# 允许 SSH (修改默认端口更安全)
sudo ufw allow 22/tcp

# 允许 HTTP/HTTPS
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

# 限制 SSH 连接速率
sudo ufw limit 22/tcp

# 启用防火墙
sudo ufw enable

# 查看状态
sudo ufw status verbose
```

### 3. DDoS 防护

```nginx
# nginx/nginx.conf
http {
    # 限制连接数
    limit_conn_zone $binary_remote_addr zone=addr:10m;
    limit_conn addr 10;
    
    # 限制请求速率
    limit_req_zone $binary_remote_addr zone=one:10m rate=10r/s;
    limit_req zone=one burst=20 nodelay;
    
    # 超时设置
    client_body_timeout 12;
    client_header_timeout 12;
    send_timeout 10;
    
    # 请求大小限制
    client_max_body_size 100M;
    client_body_buffer_size 128k;
    client_header_buffer_size 1k;
    large_client_header_buffers 4 8k;
}
```

---

## 容器安全

### 1. 非 Root 用户运行

```dockerfile
# backend/Dockerfile
FROM python:3.11-slim

# 创建非 root 用户
RUN groupadd -r ngsmodule && useradd -r -g ngsmodule ngsmodule

# 设置工作目录权限
WORKDIR /app
RUN chown -R ngsmodule:ngsmodule /app

# 切换到非 root 用户
USER ngsmodule

# 其余配置...
```

### 2. 镜像扫描

```bash
# 使用 Trivy 扫描镜像漏洞
docker run --rm -v /var/run/docker.sock:/var/run/docker.sock \
  aquasec/trivy image ngsmodule-backend:latest

# 使用 Snyk
snyk container test ngsmodule-backend:latest
```

### 3. 资源限制

```yaml
# docker-compose.yml
services:
  backend:
    deploy:
      resources:
        limits:
          cpus: '2'
          memory: 2G
          pids: 100  # 限制进程数
        reservations:
          cpus: '1'
          memory: 1G
    
    # 安全选项
    security_opt:
      - no-new-privileges:true
    
    # 只读根文件系统
    read_only: true
    tmpfs:
      - /tmp
      - /app/logs
```

### 4. 网络隔离

```yaml
# docker-compose.yml
networks:
  backend-network:
    driver: bridge
    internal: true  # 内部网络，无法访问外网
  
  frontend-network:
    driver: bridge

services:
  postgres:
    networks:
      - backend-network  # 仅后端可访问
  
  backend:
    networks:
      - backend-network
      - frontend-network
  
  frontend:
    networks:
      - frontend-network
```

---

## 密钥管理

### 1. 环境变量

```bash
# ❌ 不要在代码中硬编码密钥
# SECRET_KEY = "my-secret-key"

# ✅ 使用环境变量
# backend/app/core/config.py
from pydantic_settings import BaseSettings

class Settings(BaseSettings):
    SECRET_KEY: str
    JWT_SECRET_KEY: str
    POSTGRES_PASSWORD: str
    
    class Config:
        env_file = ".env"
        case_sensitive = True
```

### 2. Docker Secrets (生产推荐)

```yaml
# docker-compose.yml
secrets:
  postgres_password:
    file: ./secrets/postgres_password.txt
  jwt_secret:
    file: ./secrets/jwt_secret.txt

services:
  backend:
    secrets:
      - postgres_password
      - jwt_secret
    environment:
      POSTGRES_PASSWORD_FILE: /run/secrets/postgres_password
      JWT_SECRET_KEY_FILE: /run/secrets/jwt_secret
```

```python
# backend/app/core/config.py
import os

def read_secret(secret_name: str, default: str = None) -> str:
    """从 Docker secret 或环境变量读取密钥"""
    secret_file = os.getenv(f"{secret_name}_FILE")
    if secret_file and os.path.exists(secret_file):
        with open(secret_file) as f:
            return f.read().strip()
    return os.getenv(secret_name, default)

class Settings(BaseSettings):
    SECRET_KEY: str = read_secret("SECRET_KEY")
    JWT_SECRET_KEY: str = read_secret("JWT_SECRET_KEY")
    POSTGRES_PASSWORD: str = read_secret("POSTGRES_PASSWORD")
```

### 3. 密钥轮换

```python
# backend/app/core/key_rotation.py
from datetime import datetime, timedelta

class KeyRotation:
    def __init__(self):
        self.current_key = os.getenv("JWT_SECRET_KEY")
        self.previous_key = os.getenv("JWT_SECRET_KEY_OLD")
        self.rotation_date = datetime.fromisoformat(os.getenv("KEY_ROTATION_DATE"))
    
    def get_verification_keys(self) -> list[str]:
        """返回用于验证的密钥列表 (支持新旧密钥)"""
        keys = [self.current_key]
        
        # 轮换后 30 天内继续接受旧密钥
        if self.previous_key and datetime.now() - self.rotation_date < timedelta(days=30):
            keys.append(self.previous_key)
        
        return keys
    
    def get_signing_key(self) -> str:
        """返回用于签名的密钥 (始终使用当前密钥)"""
        return self.current_key
```

---

## 安全监控

### 1. 日志记录

```python
# backend/app/core/logging_config.py
import logging
from logging.handlers import RotatingFileHandler
import json

class SecurityLogger:
    def __init__(self):
        self.logger = logging.getLogger("security")
        self.logger.setLevel(logging.INFO)
        
        # 文件处理器
        handler = RotatingFileHandler(
            "logs/security.log",
            maxBytes=10485760,  # 10MB
            backupCount=10
        )
        handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(handler)
    
    def log_failed_login(self, username: str, ip: str, reason: str):
        self.logger.warning(json.dumps({
            "event": "failed_login",
            "username": username,
            "ip": ip,
            "reason": reason,
            "timestamp": datetime.now().isoformat()
        }))
    
    def log_suspicious_activity(self, user_id: str, action: str, details: dict):
        self.logger.warning(json.dumps({
            "event": "suspicious_activity",
            "user_id": user_id,
            "action": action,
            "details": details,
            "timestamp": datetime.now().isoformat()
        }))

security_logger = SecurityLogger()
```

### 2. 异常检测

```python
# backend/app/utils/anomaly_detection.py
from datetime import datetime, timedelta
from collections import defaultdict

class AnomalyDetector:
    def __init__(self, redis_client):
        self.redis = redis_client
    
    async def check_unusual_activity(self, user_id: str, action: str) -> bool:
        """检测异常活动模式"""
        key = f"activity:{user_id}:{action}"
        now = datetime.now()
        
        # 记录活动
        self.redis.zadd(key, {now.isoformat(): now.timestamp()})
        self.redis.expire(key, 3600)  # 1 小时后过期
        
        # 获取最近 1 小时的活动次数
        hour_ago = (now - timedelta(hours=1)).timestamp()
        count = self.redis.zcount(key, hour_ago, now.timestamp())
        
        # 阈值检测
        thresholds = {
            "login": 10,
            "api_call": 1000,
            "file_upload": 50,
            "delete": 20
        }
        
        if count > thresholds.get(action, 100):
            security_logger.log_suspicious_activity(
                user_id,
                action,
                {"count": count, "threshold": thresholds.get(action)}
            )
            return True
        
        return False
```

### 3. 告警集成

```python
# backend/app/utils/alerts.py
import requests

class SecurityAlerts:
    def __init__(self, webhook_url: str):
        self.webhook_url = webhook_url
    
    def send_alert(self, severity: str, title: str, message: str):
        """发送安全告警到 Slack/Teams"""
        payload = {
            "severity": severity,
            "title": title,
            "message": message,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            requests.post(self.webhook_url, json=payload, timeout=5)
        except Exception as e:
            logger.error(f"Failed to send alert: {e}")

# 使用示例
alerts = SecurityAlerts(settings.SECURITY_WEBHOOK_URL)

# 检测到可疑活动时
if await anomaly_detector.check_unusual_activity(user_id, "login"):
    alerts.send_alert(
        severity="high",
        title="Suspicious Login Activity",
        message=f"User {user_id} has attempted {count} logins in the past hour"
    )
```

---

## 安全检查清单

### 部署前检查

- [ ] **认证和授权**
  - [ ] JWT 密钥已更改 (至少 32 字符)
  - [ ] Token 过期时间配置合理
  - [ ] 实施了 Token 黑名单机制
  - [ ] 权限控制已正确配置

- [ ] **密码策略**
  - [ ] 密码强度验证已启用
  - [ ] 使用 bcrypt/argon2 哈希
  - [ ] 登录尝试限制已配置
  - [ ] 密码历史记录 (防止重复使用)

- [ ] **API 安全**
  - [ ] 速率限制已启用
  - [ ] CORS 配置仅允许可信域名
  - [ ] 输入验证已在所有端点实施
  - [ ] SQL 注入防护 (使用 ORM)

- [ ] **数据库安全**
  - [ ] 数据库密码已更改 (强密码)
  - [ ] 数据库端口未暴露到公网
  - [ ] 敏感数据已加密
  - [ ] 审计日志已启用

- [ ] **文件上传安全**
  - [ ] 文件类型验证 (扩展名 + MIME)
  - [ ] 文件大小限制
  - [ ] 文件名清理 (防止路径遍历)
  - [ ] 病毒扫描 (可选)

- [ ] **网络安全**
  - [ ] HTTPS 强制启用
  - [ ] 安全头部已配置
  - [ ] 防火墙规则已设置
  - [ ] DDoS 防护已配置

- [ ] **容器安全**
  - [ ] 所有容器以非 root 用户运行
  - [ ] 镜像漏洞已扫描
  - [ ] 资源限制已配置
  - [ ] 网络隔离已实施

- [ ] **密钥管理**
  - [ ] 无硬编码密钥
  - [ ] 环境变量或 Docker Secrets
  - [ ] 密钥轮换计划
  - [ ] .env 文件在 .gitignore 中

- [ ] **监控和日志**
  - [ ] 安全事件日志记录
  - [ ] 异常活动检测
  - [ ] 告警通知已配置
  - [ ] 日志定期审查

### 定期维护

- [ ] **每周**
  - [ ] 审查安全日志
  - [ ] 检查失败的登录尝试
  - [ ] 验证备份完整性

- [ ] **每月**
  - [ ] 更新依赖包
  - [ ] 扫描容器漏洞
  - [ ] 审查用户权限

- [ ] **每季度**
  - [ ] 轮换密钥
  - [ ] 安全审计
  - [ ] 渗透测试

---

## 应急响应

### 安全事件响应流程

1. **检测和确认**
   - 监控告警触发
   - 确认事件真实性
   - 评估影响范围

2. **遏制**
   - 隔离受影响的系统
   - 暂停可疑账户
   - 阻止攻击流量

3. **调查**
   - 收集日志和证据
   - 分析攻击向量
   - 确定数据泄露范围

4. **恢复**
   - 修复漏洞
   - 恢复受影响的服务
   - 重置受影响的密钥

5. **后续**
   - 更新安全策略
   - 加强监控
   - 团队培训

### 紧急联系方式

```bash
# 立即撤销所有 Token
docker-compose exec backend python scripts/revoke_all_tokens.py

# 禁用所有用户登录
docker-compose exec backend python scripts/emergency_lockdown.py

# 查看最近的安全日志
docker-compose exec backend tail -n 100 /app/logs/security.log

# 导出审计日志
docker-compose exec postgres psql -U ngsmodule -d ngsmodule -c \
  "COPY (SELECT * FROM audit_logs WHERE created_at > NOW() - INTERVAL '24 hours') TO STDOUT CSV HEADER" > audit_export.csv
```

---

**最后更新**: 2025-01-22  
**版本**: 1.0.0

**重要提示**: 安全是一个持续的过程，而非一次性任务。定期审查和更新安全措施至关重要。
