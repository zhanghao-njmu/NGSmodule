"""
Observability: Prometheus metrics and Sentry error tracking.

Both integrations are optional and gracefully degrade when their
dependencies are missing.
"""

import logging
from typing import Optional

from fastapi import FastAPI

from app.core.config import settings

logger = logging.getLogger(__name__)


# ----------------------------------------------------------------------
# Sentry
# ----------------------------------------------------------------------


def init_sentry() -> bool:
    """Initialize Sentry SDK if SENTRY_DSN is configured.

    Returns True if initialized.
    """
    dsn = getattr(settings, "SENTRY_DSN", None)
    if not dsn:
        return False

    try:
        import sentry_sdk
        from sentry_sdk.integrations.celery import CeleryIntegration
        from sentry_sdk.integrations.fastapi import FastApiIntegration
        from sentry_sdk.integrations.sqlalchemy import SqlalchemyIntegration
        from sentry_sdk.integrations.starlette import StarletteIntegration
    except ImportError:
        logger.warning(
            "SENTRY_DSN configured but sentry-sdk is not installed. " "Run: pip install 'sentry-sdk[fastapi]'"
        )
        return False

    sample_rate = getattr(settings, "SENTRY_TRACES_SAMPLE_RATE", 0.1)
    environment = getattr(settings, "SENTRY_ENVIRONMENT", "production")

    sentry_sdk.init(
        dsn=dsn,
        traces_sample_rate=sample_rate,
        environment=environment,
        release=settings.APP_VERSION,
        integrations=[
            FastApiIntegration(),
            StarletteIntegration(),
            SqlalchemyIntegration(),
            CeleryIntegration(),
        ],
        # Don't send PII by default; let operators opt-in via env
        send_default_pii=getattr(settings, "SENTRY_SEND_PII", False),
    )

    logger.info(f"Sentry initialized for environment={environment}")
    return True


# ----------------------------------------------------------------------
# Prometheus
# ----------------------------------------------------------------------


def init_prometheus(app: FastAPI) -> bool:
    """Mount Prometheus metrics endpoint at /metrics.

    Uses prometheus-fastapi-instrumentator for automatic HTTP metrics
    collection (request count, latency histogram, in-progress, etc.).

    Returns True if mounted.
    """
    if not getattr(settings, "PROMETHEUS_ENABLED", True):
        return False

    try:
        from prometheus_fastapi_instrumentator import Instrumentator
    except ImportError:
        logger.info(
            "Prometheus instrumentation skipped: prometheus-fastapi-instrumentator "
            "is not installed. Run: pip install prometheus-fastapi-instrumentator"
        )
        return False

    instrumentator = Instrumentator(
        should_group_status_codes=True,
        should_ignore_untemplated=True,
        should_respect_env_var=True,
        env_var_name="ENABLE_METRICS",
        excluded_handlers=["/metrics", "/health"],
    )

    instrumentator.instrument(app).expose(
        app,
        endpoint="/metrics",
        include_in_schema=False,
        tags=["Observability"],
    )

    logger.info("Prometheus metrics endpoint mounted at /metrics")
    return True


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def init_observability(app: Optional[FastAPI] = None) -> dict:
    """Initialize all observability backends.

    Returns a dict with the status of each integration.
    """
    status = {
        "sentry": init_sentry(),
        "prometheus": False,
    }
    if app is not None:
        status["prometheus"] = init_prometheus(app)
    return status
