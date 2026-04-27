"""
Application Configuration
"""
from typing import List, Optional
from pydantic_settings import BaseSettings
from pydantic import AnyHttpUrl, validator


class Settings(BaseSettings):
    """Application settings"""

    # Application
    APP_NAME: str = "NGSmodule"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False
    SECRET_KEY: str

    # API
    API_V1_PREFIX: str = "/api/v1"
    BACKEND_CORS_ORIGINS: List[AnyHttpUrl] = []

    @validator("BACKEND_CORS_ORIGINS", pre=True)
    def assemble_cors_origins(cls, v: str | List[str]) -> List[str] | str:
        if isinstance(v, str) and not v.startswith("["):
            return [i.strip() for i in v.split(",")]
        elif isinstance(v, (list, str)):
            return v
        raise ValueError(v)

    # Database
    DATABASE_URL: str
    DATABASE_POOL_SIZE: int = 20
    DATABASE_MAX_OVERFLOW: int = 0

    # Redis
    REDIS_URL: str = "redis://localhost:6379/0"

    # Celery
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/0"

    # MinIO / S3
    MINIO_URL: str = "localhost:9000"
    MINIO_ACCESS_KEY: str
    MINIO_SECRET_KEY: str
    MINIO_BUCKET: str = "ngsmodule"
    MINIO_SECURE: bool = False

    # JWT
    JWT_SECRET_KEY: str
    JWT_ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    REFRESH_TOKEN_EXPIRE_DAYS: int = 7

    # File Storage
    UPLOAD_DIR: str = "/data/uploads"
    MAX_UPLOAD_SIZE: int = 53687091200  # 50GB
    ALLOWED_EXTENSIONS: List[str] = [
        ".fastq", ".fastq.gz", ".fq", ".fq.gz",
        ".bam", ".sam", ".vcf", ".vcf.gz"
    ]

    # NGS Pipeline
    NGS_PIPELINE_DIR: str
    NGS_WORK_DIR: str = "/data/ngsmodule_work"
    NGS_RESULT_DIR: str = "/data/ngsmodule_results"

    # Logging
    LOG_LEVEL: str = "INFO"
    LOG_FILE: Optional[str] = None

    # Backup
    BACKUP_DIR: str = "/data/backups"
    BACKUP_RETENTION_DAYS: int = 30

    # Alerts (thresholds)
    ALERT_DISK_WARNING_PERCENT: float = 80.0
    ALERT_DISK_CRITICAL_PERCENT: float = 90.0
    ALERT_MEMORY_WARNING_PERCENT: float = 85.0
    ALERT_MEMORY_CRITICAL_PERCENT: float = 95.0

    # Observability
    SENTRY_DSN: Optional[str] = None
    SENTRY_ENVIRONMENT: str = "production"
    SENTRY_TRACES_SAMPLE_RATE: float = 0.1
    SENTRY_SEND_PII: bool = False
    PROMETHEUS_ENABLED: bool = True

    # Rate limiting
    RATE_LIMIT_STORAGE: Optional[str] = None  # Override REDIS_URL for slowapi

    # AI provider
    AI_PROVIDER: str = "mock"  # mock | claude | openai
    ANTHROPIC_API_KEY: Optional[str] = None
    ANTHROPIC_MODEL: str = "claude-sonnet-4-6"
    OPENAI_API_KEY: Optional[str] = None
    OPENAI_MODEL: str = "gpt-4o-mini"

    class Config:
        env_file = ".env"
        case_sensitive = True


settings = Settings()
