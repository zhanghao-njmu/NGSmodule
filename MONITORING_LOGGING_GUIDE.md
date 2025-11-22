# NGSmodule 监控和日志配置指南

本文档提供 NGSmodule 生产环境的监控和日志配置方案。

## 📋 目录

- [架构概览](#架构概览)
- [Prometheus 指标收集](#prometheus-指标收集)
- [Grafana 可视化](#grafana-可视化)
- [日志管理](#日志管理)
- [告警配置](#告警配置)
- [健康检查](#健康检查)
- [性能监控](#性能监控)
- [故障排查](#故障排查)

---

## 架构概览

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│  Backend    │─────▶│ Prometheus  │─────▶│  Grafana    │
│  Frontend   │      │  (Metrics)  │      │ (Dashboard) │
│  Database   │      └─────────────┘      └─────────────┘
└─────────────┘              │
       │                     │
       │              ┌──────▼──────┐
       │              │ AlertManager│
       │              └─────────────┘
       │
       ▼
┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   Logs      │─────▶│  Loki/ELK   │─────▶│  Grafana    │
│  (stdout)   │      │   (Storage) │      │  (Viewer)   │
└─────────────┘      └─────────────┘      └─────────────┘
```

---

## Prometheus 指标收集

### 1. Backend 指标导出

安装 Prometheus 客户端:

```bash
# backend/requirements.txt
prometheus-client==0.19.0
prometheus-fastapi-instrumentator==6.1.0
```

配置指标端点:

```python
# backend/app/main.py
from prometheus_fastapi_instrumentator import Instrumentator
from prometheus_client import Counter, Histogram, Gauge
import time

# 初始化 Instrumentator
instrumentator = Instrumentator(
    should_group_status_codes=True,
    should_ignore_untemplated=True,
    should_respect_env_var=True,
    should_instrument_requests_inprogress=True,
    excluded_handlers=["/metrics", "/health"],
    env_var_name="ENABLE_METRICS",
    inprogress_name="http_requests_inprogress",
    inprogress_labels=True,
)

# 自定义指标
request_duration = Histogram(
    'http_request_duration_seconds',
    'HTTP request duration in seconds',
    ['method', 'endpoint', 'status']
)

db_query_duration = Histogram(
    'db_query_duration_seconds',
    'Database query duration in seconds',
    ['query_type']
)

active_users = Gauge(
    'active_users_total',
    'Number of active users'
)

file_uploads = Counter(
    'file_uploads_total',
    'Total number of file uploads',
    ['status']
)

pipeline_tasks = Gauge(
    'pipeline_tasks_total',
    'Number of pipeline tasks',
    ['status']
)

# 中间件
@app.middleware("http")
async def prometheus_middleware(request: Request, call_next):
    start_time = time.time()
    response = await call_next(request)
    duration = time.time() - start_time
    
    request_duration.labels(
        method=request.method,
        endpoint=request.url.path,
        status=response.status_code
    ).observe(duration)
    
    return response

# 启用指标
instrumentator.instrument(app).expose(app, endpoint="/metrics")

# 自定义指标更新
@app.on_event("startup")
async def update_metrics():
    """定期更新自定义指标"""
    import asyncio
    
    async def update_loop():
        while True:
            try:
                db = SessionLocal()
                
                # 活跃用户数
                active_count = db.query(User).filter(
                    User.is_active == True,
                    User.last_login > datetime.now() - timedelta(days=30)
                ).count()
                active_users.set(active_count)
                
                # 任务状态统计
                task_stats = db.query(
                    PipelineTask.status,
                    func.count(PipelineTask.id)
                ).group_by(PipelineTask.status).all()
                
                for status, count in task_stats:
                    pipeline_tasks.labels(status=status).set(count)
                
                db.close()
            except Exception as e:
                logger.error(f"Error updating metrics: {e}")
            
            await asyncio.sleep(60)  # 每分钟更新
    
    asyncio.create_task(update_loop())
```

### 2. PostgreSQL 导出器

```yaml
# docker-compose.yml
services:
  postgres-exporter:
    image: prometheuscommunity/postgres-exporter:latest
    environment:
      DATA_SOURCE_NAME: "postgresql://${POSTGRES_USER}:${POSTGRES_PASSWORD}@postgres:5432/${POSTGRES_DB}?sslmode=disable"
    ports:
      - "9187:9187"
    networks:
      - backend-network
    depends_on:
      - postgres
```

### 3. Redis 导出器

```yaml
# docker-compose.yml
services:
  redis-exporter:
    image: oliver006/redis_exporter:latest
    environment:
      REDIS_ADDR: redis:6379
      REDIS_PASSWORD: ${REDIS_PASSWORD}
    ports:
      - "9121:9121"
    networks:
      - backend-network
    depends_on:
      - redis
```

### 4. Nginx 导出器

```yaml
# docker-compose.yml
services:
  nginx-exporter:
    image: nginx/nginx-prometheus-exporter:latest
    command:
      - '-nginx.scrape-uri=http://nginx:80/stub_status'
    ports:
      - "9113:9113"
    networks:
      - frontend-network
    depends_on:
      - nginx
```

配置 Nginx stub_status:

```nginx
# frontend/nginx.conf
server {
    listen 80;
    
    # Metrics endpoint
    location /stub_status {
        stub_status on;
        access_log off;
        allow 172.16.0.0/12;  # Docker 内部网络
        deny all;
    }
}
```

### 5. Prometheus 配置

```yaml
# monitoring/prometheus/prometheus.yml
global:
  scrape_interval: 15s
  evaluation_interval: 15s
  external_labels:
    cluster: 'ngsmodule-production'
    environment: 'production'

# 告警规则
rule_files:
  - '/etc/prometheus/rules/*.yml'

# 抓取配置
scrape_configs:
  # Backend API
  - job_name: 'backend'
    static_configs:
      - targets: ['backend:8000']
    metrics_path: '/metrics'
    scrape_interval: 10s

  # PostgreSQL
  - job_name: 'postgres'
    static_configs:
      - targets: ['postgres-exporter:9187']
    scrape_interval: 30s

  # Redis
  - job_name: 'redis'
    static_configs:
      - targets: ['redis-exporter:9121']
    scrape_interval: 30s

  # Nginx
  - job_name: 'nginx'
    static_configs:
      - targets: ['nginx-exporter:9113']
    scrape_interval: 30s

  # Node Exporter (系统指标)
  - job_name: 'node'
    static_configs:
      - targets: ['node-exporter:9100']
    scrape_interval: 30s

  # Celery
  - job_name: 'celery'
    static_configs:
      - targets: ['celery-exporter:9808']
    scrape_interval: 30s
```

Docker Compose 配置:

```yaml
# docker-compose.yml
services:
  prometheus:
    image: prom/prometheus:latest
    volumes:
      - ./monitoring/prometheus/prometheus.yml:/etc/prometheus/prometheus.yml:ro
      - ./monitoring/prometheus/rules:/etc/prometheus/rules:ro
      - prometheus-data:/prometheus
    command:
      - '--config.file=/etc/prometheus/prometheus.yml'
      - '--storage.tsdb.path=/prometheus'
      - '--storage.tsdb.retention.time=30d'
      - '--web.console.libraries=/usr/share/prometheus/console_libraries'
      - '--web.console.templates=/usr/share/prometheus/consoles'
    ports:
      - "9090:9090"
    networks:
      - backend-network
    restart: unless-stopped

  node-exporter:
    image: prom/node-exporter:latest
    command:
      - '--path.procfs=/host/proc'
      - '--path.sysfs=/host/sys'
      - '--collector.filesystem.mount-points-exclude=^/(sys|proc|dev|host|etc)($$|/)'
    volumes:
      - /proc:/host/proc:ro
      - /sys:/host/sys:ro
      - /:/rootfs:ro
    ports:
      - "9100:9100"
    networks:
      - backend-network

volumes:
  prometheus-data:
```

---

## Grafana 可视化

### 1. Grafana 配置

```yaml
# docker-compose.yml
services:
  grafana:
    image: grafana/grafana:latest
    environment:
      - GF_SECURITY_ADMIN_USER=admin
      - GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_PASSWORD}
      - GF_SERVER_ROOT_URL=https://grafana.yourdomain.com
      - GF_INSTALL_PLUGINS=grafana-clock-panel,grafana-simple-json-datasource
    volumes:
      - grafana-data:/var/lib/grafana
      - ./monitoring/grafana/provisioning:/etc/grafana/provisioning:ro
      - ./monitoring/grafana/dashboards:/var/lib/grafana/dashboards:ro
    ports:
      - "3001:3000"
    networks:
      - backend-network
    depends_on:
      - prometheus
    restart: unless-stopped

volumes:
  grafana-data:
```

### 2. 数据源配置

```yaml
# monitoring/grafana/provisioning/datasources/prometheus.yml
apiVersion: 1

datasources:
  - name: Prometheus
    type: prometheus
    access: proxy
    url: http://prometheus:9090
    isDefault: true
    editable: false
```

### 3. 仪表板配置

```yaml
# monitoring/grafana/provisioning/dashboards/dashboards.yml
apiVersion: 1

providers:
  - name: 'NGSmodule Dashboards'
    orgId: 1
    folder: ''
    type: file
    disableDeletion: false
    updateIntervalSeconds: 10
    allowUiUpdates: true
    options:
      path: /var/lib/grafana/dashboards
```

### 4. 示例仪表板

创建 `monitoring/grafana/dashboards/ngsmodule-overview.json`:

```json
{
  "dashboard": {
    "title": "NGSmodule Overview",
    "panels": [
      {
        "title": "Request Rate",
        "targets": [
          {
            "expr": "rate(http_requests_total[5m])",
            "legendFormat": "{{method}} {{endpoint}}"
          }
        ]
      },
      {
        "title": "Response Time",
        "targets": [
          {
            "expr": "histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m]))",
            "legendFormat": "p95"
          },
          {
            "expr": "histogram_quantile(0.99, rate(http_request_duration_seconds_bucket[5m]))",
            "legendFormat": "p99"
          }
        ]
      },
      {
        "title": "Error Rate",
        "targets": [
          {
            "expr": "rate(http_requests_total{status=~\"5..\"}[5m])",
            "legendFormat": "5xx errors"
          }
        ]
      },
      {
        "title": "Database Query Duration",
        "targets": [
          {
            "expr": "rate(db_query_duration_seconds_sum[5m]) / rate(db_query_duration_seconds_count[5m])",
            "legendFormat": "{{query_type}}"
          }
        ]
      },
      {
        "title": "Active Users",
        "targets": [
          {
            "expr": "active_users_total"
          }
        ]
      },
      {
        "title": "Pipeline Tasks by Status",
        "targets": [
          {
            "expr": "pipeline_tasks_total",
            "legendFormat": "{{status}}"
          }
        ]
      },
      {
        "title": "CPU Usage",
        "targets": [
          {
            "expr": "100 - (avg by (instance) (irate(node_cpu_seconds_total{mode=\"idle\"}[5m])) * 100)"
          }
        ]
      },
      {
        "title": "Memory Usage",
        "targets": [
          {
            "expr": "(node_memory_MemTotal_bytes - node_memory_MemAvailable_bytes) / node_memory_MemTotal_bytes * 100"
          }
        ]
      }
    ]
  }
}
```

---

## 日志管理

### 方案 1: Grafana Loki (推荐)

#### 1. Loki 配置

```yaml
# docker-compose.yml
services:
  loki:
    image: grafana/loki:latest
    ports:
      - "3100:3100"
    volumes:
      - ./monitoring/loki/loki-config.yml:/etc/loki/local-config.yaml:ro
      - loki-data:/loki
    command: -config.file=/etc/loki/local-config.yaml
    networks:
      - backend-network
    restart: unless-stopped

  promtail:
    image: grafana/promtail:latest
    volumes:
      - ./monitoring/promtail/promtail-config.yml:/etc/promtail/config.yml:ro
      - /var/lib/docker/containers:/var/lib/docker/containers:ro
      - /var/log:/var/log:ro
    command: -config.file=/etc/promtail/config.yml
    networks:
      - backend-network
    depends_on:
      - loki
    restart: unless-stopped

volumes:
  loki-data:
```

#### 2. Loki 配置文件

```yaml
# monitoring/loki/loki-config.yml
auth_enabled: false

server:
  http_listen_port: 3100

ingester:
  lifecycler:
    address: 127.0.0.1
    ring:
      kvstore:
        store: inmemory
      replication_factor: 1
  chunk_idle_period: 5m
  chunk_retain_period: 30s

schema_config:
  configs:
    - from: 2024-01-01
      store: boltdb-shipper
      object_store: filesystem
      schema: v11
      index:
        prefix: index_
        period: 24h

storage_config:
  boltdb_shipper:
    active_index_directory: /loki/index
    cache_location: /loki/cache
    shared_store: filesystem
  filesystem:
    directory: /loki/chunks

limits_config:
  enforce_metric_name: false
  reject_old_samples: true
  reject_old_samples_max_age: 168h
  retention_period: 720h  # 30 天

chunk_store_config:
  max_look_back_period: 0s

table_manager:
  retention_deletes_enabled: true
  retention_period: 720h
```

#### 3. Promtail 配置

```yaml
# monitoring/promtail/promtail-config.yml
server:
  http_listen_port: 9080
  grpc_listen_port: 0

positions:
  filename: /tmp/positions.yaml

clients:
  - url: http://loki:3100/loki/api/v1/push

scrape_configs:
  # Docker 容器日志
  - job_name: docker
    docker_sd_configs:
      - host: unix:///var/run/docker.sock
        refresh_interval: 5s
    relabel_configs:
      - source_labels: ['__meta_docker_container_name']
        regex: '/(.*)'
        target_label: 'container'
      - source_labels: ['__meta_docker_container_log_stream']
        target_label: 'logstream'
      - source_labels: ['__meta_docker_container_label_com_docker_compose_service']
        target_label: 'service'

  # 系统日志
  - job_name: system
    static_configs:
      - targets:
          - localhost
        labels:
          job: varlogs
          __path__: /var/log/*.log
```

#### 4. Grafana 添加 Loki 数据源

```yaml
# monitoring/grafana/provisioning/datasources/loki.yml
apiVersion: 1

datasources:
  - name: Loki
    type: loki
    access: proxy
    url: http://loki:3100
    editable: false
```

### 方案 2: ELK Stack

```yaml
# docker-compose.yml
services:
  elasticsearch:
    image: docker.elastic.co/elasticsearch/elasticsearch:8.11.0
    environment:
      - discovery.type=single-node
      - "ES_JAVA_OPTS=-Xms512m -Xmx512m"
      - xpack.security.enabled=false
    volumes:
      - elasticsearch-data:/usr/share/elasticsearch/data
    ports:
      - "9200:9200"
    networks:
      - backend-network

  logstash:
    image: docker.elastic.co/logstash/logstash:8.11.0
    volumes:
      - ./monitoring/logstash/logstash.conf:/usr/share/logstash/pipeline/logstash.conf:ro
    ports:
      - "5000:5000"
    networks:
      - backend-network
    depends_on:
      - elasticsearch

  kibana:
    image: docker.elastic.co/kibana/kibana:8.11.0
    environment:
      - ELASTICSEARCH_HOSTS=http://elasticsearch:9200
    ports:
      - "5601:5601"
    networks:
      - backend-network
    depends_on:
      - elasticsearch

volumes:
  elasticsearch-data:
```

### 应用日志配置

#### Backend 结构化日志

```python
# backend/app/core/logging_config.py
import logging
import json
from datetime import datetime
import sys

class JSONFormatter(logging.Formatter):
    """JSON 格式日志"""
    def format(self, record):
        log_data = {
            "timestamp": datetime.utcnow().isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
            "module": record.module,
            "function": record.funcName,
            "line": record.lineno,
        }
        
        # 添加异常信息
        if record.exc_info:
            log_data["exception"] = self.formatException(record.exc_info)
        
        # 添加额外字段
        if hasattr(record, 'user_id'):
            log_data["user_id"] = record.user_id
        if hasattr(record, 'request_id'):
            log_data["request_id"] = record.request_id
        
        return json.dumps(log_data)

def setup_logging():
    """配置日志系统"""
    # 根日志器
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    # 控制台处理器
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(JSONFormatter())
    root_logger.addHandler(console_handler)
    
    # 文件处理器 (可选)
    if os.path.exists('/app/logs'):
        from logging.handlers import RotatingFileHandler
        file_handler = RotatingFileHandler(
            '/app/logs/app.log',
            maxBytes=10485760,  # 10MB
            backupCount=10
        )
        file_handler.setFormatter(JSONFormatter())
        root_logger.addHandler(file_handler)
    
    return root_logger

# backend/app/main.py
from app.core.logging_config import setup_logging

logger = setup_logging()
```

#### 请求日志中间件

```python
# backend/app/middleware/logging.py
import uuid
import time
import logging
from fastapi import Request

logger = logging.getLogger(__name__)

@app.middleware("http")
async def logging_middleware(request: Request, call_next):
    # 生成请求 ID
    request_id = str(uuid.uuid4())
    request.state.request_id = request_id
    
    # 记录请求
    start_time = time.time()
    logger.info(
        "Request started",
        extra={
            "request_id": request_id,
            "method": request.method,
            "path": request.url.path,
            "client_ip": request.client.host,
            "user_agent": request.headers.get("user-agent")
        }
    )
    
    # 处理请求
    try:
        response = await call_next(request)
        duration = time.time() - start_time
        
        # 记录响应
        logger.info(
            "Request completed",
            extra={
                "request_id": request_id,
                "status_code": response.status_code,
                "duration": f"{duration:.3f}s"
            }
        )
        
        # 添加请求 ID 到响应头
        response.headers["X-Request-ID"] = request_id
        return response
        
    except Exception as e:
        duration = time.time() - start_time
        logger.error(
            "Request failed",
            extra={
                "request_id": request_id,
                "duration": f"{duration:.3f}s",
                "error": str(e)
            },
            exc_info=True
        )
        raise
```

---

## 告警配置

### 1. Alertmanager 配置

```yaml
# docker-compose.yml
services:
  alertmanager:
    image: prom/alertmanager:latest
    volumes:
      - ./monitoring/alertmanager/alertmanager.yml:/etc/alertmanager/alertmanager.yml:ro
    command:
      - '--config.file=/etc/alertmanager/alertmanager.yml'
    ports:
      - "9093:9093"
    networks:
      - backend-network
    restart: unless-stopped
```

```yaml
# monitoring/alertmanager/alertmanager.yml
global:
  resolve_timeout: 5m
  smtp_smarthost: 'smtp.gmail.com:587'
  smtp_from: 'alerts@yourdomain.com'
  smtp_auth_username: 'alerts@yourdomain.com'
  smtp_auth_password: 'your_password'

route:
  group_by: ['alertname', 'cluster']
  group_wait: 10s
  group_interval: 10s
  repeat_interval: 12h
  receiver: 'default'
  routes:
    - match:
        severity: critical
      receiver: 'critical'
    - match:
        severity: warning
      receiver: 'warning'

receivers:
  - name: 'default'
    email_configs:
      - to: 'team@yourdomain.com'

  - name: 'critical'
    email_configs:
      - to: 'oncall@yourdomain.com'
        headers:
          Subject: '[CRITICAL] {{ .GroupLabels.alertname }}'
    slack_configs:
      - api_url: 'YOUR_SLACK_WEBHOOK_URL'
        channel: '#alerts-critical'
        title: 'Critical Alert: {{ .GroupLabels.alertname }}'
        text: '{{ range .Alerts }}{{ .Annotations.description }}{{ end }}'

  - name: 'warning'
    email_configs:
      - to: 'team@yourdomain.com'
    slack_configs:
      - api_url: 'YOUR_SLACK_WEBHOOK_URL'
        channel: '#alerts-warning'
```

### 2. 告警规则

```yaml
# monitoring/prometheus/rules/alerts.yml
groups:
  - name: ngsmodule_alerts
    interval: 30s
    rules:
      # 高错误率
      - alert: HighErrorRate
        expr: rate(http_requests_total{status=~"5.."}[5m]) > 0.05
        for: 5m
        labels:
          severity: critical
        annotations:
          summary: "High error rate detected"
          description: "Error rate is {{ $value | humanizePercentage }} for {{ $labels.instance }}"

      # 响应时间过长
      - alert: HighResponseTime
        expr: histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m])) > 1
        for: 10m
        labels:
          severity: warning
        annotations:
          summary: "High response time detected"
          description: "95th percentile response time is {{ $value }}s"

      # 数据库连接数过高
      - alert: HighDatabaseConnections
        expr: pg_stat_database_numbackends > 80
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High number of database connections"
          description: "{{ $value }} active connections to {{ $labels.datname }}"

      # 磁盘空间不足
      - alert: LowDiskSpace
        expr: (node_filesystem_avail_bytes / node_filesystem_size_bytes) * 100 < 15
        for: 5m
        labels:
          severity: critical
        annotations:
          summary: "Low disk space on {{ $labels.instance }}"
          description: "Only {{ $value | humanize }}% disk space left"

      # 内存使用率过高
      - alert: HighMemoryUsage
        expr: (1 - (node_memory_MemAvailable_bytes / node_memory_MemTotal_bytes)) * 100 > 90
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "High memory usage on {{ $labels.instance }}"
          description: "Memory usage is {{ $value | humanize }}%"

      # CPU 使用率过高
      - alert: HighCPUUsage
        expr: 100 - (avg by (instance) (irate(node_cpu_seconds_total{mode="idle"}[5m])) * 100) > 80
        for: 10m
        labels:
          severity: warning
        annotations:
          summary: "High CPU usage on {{ $labels.instance }}"
          description: "CPU usage is {{ $value | humanize }}%"

      # Redis 内存使用率过高
      - alert: RedisHighMemory
        expr: (redis_memory_used_bytes / redis_memory_max_bytes) * 100 > 80
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "Redis memory usage high"
          description: "Redis is using {{ $value | humanize }}% of max memory"

      # Celery 队列堆积
      - alert: CeleryQueueBacklog
        expr: celery_queue_length > 100
        for: 10m
        labels:
          severity: warning
        annotations:
          summary: "Celery queue backlog"
          description: "{{ $value }} tasks waiting in queue {{ $labels.queue_name }}"

      # 服务不可用
      - alert: ServiceDown
        expr: up == 0
        for: 1m
        labels:
          severity: critical
        annotations:
          summary: "Service {{ $labels.job }} is down"
          description: "{{ $labels.instance }} has been down for more than 1 minute"
```

---

## 健康检查

### 1. Backend 健康检查端点

```python
# backend/app/api/v1/health.py
from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from app.db.session import get_db
from app.core.redis import redis_client
import psutil

router = APIRouter()

@router.get("/health")
async def health_check():
    """基本健康检查"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat()
    }

@router.get("/health/detailed")
async def detailed_health_check(db: Session = Depends(get_db)):
    """详细健康检查"""
    checks = {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "checks": {}
    }
    
    # 数据库检查
    try:
        db.execute("SELECT 1")
        checks["checks"]["database"] = {"status": "healthy"}
    except Exception as e:
        checks["status"] = "unhealthy"
        checks["checks"]["database"] = {
            "status": "unhealthy",
            "error": str(e)
        }
    
    # Redis 检查
    try:
        redis_client.ping()
        checks["checks"]["redis"] = {"status": "healthy"}
    except Exception as e:
        checks["status"] = "unhealthy"
        checks["checks"]["redis"] = {
            "status": "unhealthy",
            "error": str(e)
        }
    
    # 磁盘空间检查
    disk_usage = psutil.disk_usage('/')
    checks["checks"]["disk"] = {
        "status": "healthy" if disk_usage.percent < 85 else "warning",
        "usage_percent": disk_usage.percent,
        "free_gb": disk_usage.free / (1024**3)
    }
    
    # 内存检查
    memory = psutil.virtual_memory()
    checks["checks"]["memory"] = {
        "status": "healthy" if memory.percent < 85 else "warning",
        "usage_percent": memory.percent,
        "available_gb": memory.available / (1024**3)
    }
    
    return checks

@router.get("/health/ready")
async def readiness_check(db: Session = Depends(get_db)):
    """Kubernetes readiness probe"""
    try:
        db.execute("SELECT 1")
        redis_client.ping()
        return {"ready": True}
    except Exception as e:
        return {"ready": False, "error": str(e)}, 503

@router.get("/health/live")
async def liveness_check():
    """Kubernetes liveness probe"""
    return {"alive": True}
```

---

## 性能监控

### 1. 应用性能监控 (APM)

使用 OpenTelemetry:

```python
# backend/requirements.txt
opentelemetry-api==1.21.0
opentelemetry-sdk==1.21.0
opentelemetry-instrumentation-fastapi==0.42b0
opentelemetry-instrumentation-sqlalchemy==0.42b0
opentelemetry-exporter-otlp==1.21.0

# backend/app/core/tracing.py
from opentelemetry import trace
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor
from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter
from opentelemetry.instrumentation.fastapi import FastAPIInstrumentor
from opentelemetry.instrumentation.sqlalchemy import SQLAlchemyInstrumentor

def setup_tracing(app):
    """配置 OpenTelemetry 追踪"""
    # 设置追踪提供者
    trace.set_tracer_provider(TracerProvider())
    tracer_provider = trace.get_tracer_provider()
    
    # OTLP 导出器
    otlp_exporter = OTLPSpanExporter(
        endpoint="http://tempo:4317",
        insecure=True
    )
    
    # 批处理 span 处理器
    span_processor = BatchSpanProcessor(otlp_exporter)
    tracer_provider.add_span_processor(span_processor)
    
    # 自动注入 FastAPI
    FastAPIInstrumentor.instrument_app(app)
    
    # 自动注入 SQLAlchemy
    from app.db.session import engine
    SQLAlchemyInstrumentor().instrument(engine=engine)

# backend/app/main.py
from app.core.tracing import setup_tracing

setup_tracing(app)
```

---

## 故障排查

### 常用查询

#### Prometheus 查询

```promql
# 每秒请求数
rate(http_requests_total[5m])

# 95th 百分位响应时间
histogram_quantile(0.95, rate(http_request_duration_seconds_bucket[5m]))

# 错误率
rate(http_requests_total{status=~"5.."}[5m]) / rate(http_requests_total[5m])

# 数据库查询延迟
rate(db_query_duration_seconds_sum[5m]) / rate(db_query_duration_seconds_count[5m])

# 内存使用率
(node_memory_MemTotal_bytes - node_memory_MemAvailable_bytes) / node_memory_MemTotal_bytes * 100
```

#### Loki 查询

```logql
# 查找错误日志
{container="backend"} |= "ERROR"

# 统计错误数
sum(rate({container="backend"} |= "ERROR" [5m]))

# 查找特定用户的请求
{container="backend"} | json | user_id="123"

# 查找慢查询
{container="backend"} | json | duration > 1s
```

---

**最后更新**: 2025-01-22  
**版本**: 1.0.0
