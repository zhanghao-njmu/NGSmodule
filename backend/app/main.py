"""
Main application entry point
"""
import logging
from contextlib import asynccontextmanager
from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from sqlalchemy.exc import IntegrityError, DBAPIError
from pydantic import ValidationError
from slowapi import _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from app.core.config import settings
from app.core.database import init_db
from app.core.rate_limit import limiter, rate_limit_exceeded_handler
from app.api.v1 import auth, users, projects, samples, files, tasks, websocket, pipelines, results, stats, notifications, analytics, admin, ai, data_downloads, vendor_credentials

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan: startup + shutdown."""
    # Startup
    init_db()
    logger.info(f"{settings.APP_NAME} v{settings.APP_VERSION} started successfully")
    logger.info(f"API documentation: http://localhost:8000{settings.API_V1_PREFIX}/docs")

    # Initialize observability (Sentry, Prometheus) lazily so missing deps
    # don't break startup
    try:
        from app.core.observability import init_observability
        init_observability(app)
    except Exception as e:
        logger.warning(f"Observability initialization skipped: {e}")

    # Start realtime subscriber for WebSocket fanout
    realtime_sub = None
    try:
        from app.api.v1.websocket import start_realtime_subscriber
        realtime_sub = await start_realtime_subscriber()
    except Exception as e:
        logger.warning(f"Realtime subscriber not started: {e}")

    yield

    # Shutdown
    if realtime_sub is not None:
        try:
            await realtime_sub.stop()
        except Exception:
            pass
    logger.info(f"{settings.APP_NAME} shutting down...")


# Create FastAPI application
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="NGSmodule API - Enterprise-grade bioinformatics workstation",
    docs_url=f"{settings.API_V1_PREFIX}/docs",
    redoc_url=f"{settings.API_V1_PREFIX}/redoc",
    openapi_url=f"{settings.API_V1_PREFIX}/openapi.json",
    lifespan=lifespan,
)

# Add rate limiter state to app
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, rate_limit_exceeded_handler)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.BACKEND_CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Global Exception Handlers
@app.exception_handler(IntegrityError)
async def integrity_error_handler(request: Request, exc: IntegrityError):
    """Handle database integrity constraint violations"""
    logger.error(f"Database integrity error: {str(exc)}")
    return JSONResponse(
        status_code=status.HTTP_409_CONFLICT,
        content={
            "error": "Database constraint violation",
            "detail": "A unique constraint was violated. The resource may already exist.",
            "code": "INTEGRITY_ERROR"
        }
    )


@app.exception_handler(DBAPIError)
async def database_error_handler(request: Request, exc: DBAPIError):
    """Handle database connection and query errors"""
    logger.error(f"Database error: {str(exc)}", exc_info=True)
    return JSONResponse(
        status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
        content={
            "error": "Database error",
            "detail": "Unable to connect to database or execute query. Please try again later.",
            "code": "DATABASE_ERROR"
        }
    )


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle FastAPI request validation errors"""
    logger.warning(f"Validation error: {exc.errors()}")
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
        content={
            "error": "Validation error",
            "detail": exc.errors(),
            "code": "VALIDATION_ERROR"
        }
    )


@app.exception_handler(ValidationError)
async def pydantic_validation_error_handler(request: Request, exc: ValidationError):
    """Handle Pydantic validation errors"""
    logger.warning(f"Pydantic validation error: {exc.errors()}")
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
        content={
            "error": "Validation error",
            "detail": exc.errors(),
            "code": "VALIDATION_ERROR"
        }
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Catch all unhandled exceptions"""
    logger.error(f"Unhandled exception: {str(exc)}", exc_info=True)

    # Don't expose internal error details in production
    detail = str(exc) if settings.DEBUG else "An unexpected error occurred"

    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            "error": "Internal server error",
            "detail": detail,
            "code": "INTERNAL_ERROR"
        }
    )


@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "running",
        "docs": f"{settings.API_V1_PREFIX}/docs"
    }


@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy"}


# Include routers
app.include_router(
    auth.router,
    prefix=f"{settings.API_V1_PREFIX}/auth",
    tags=["Authentication"]
)

app.include_router(
    users.router,
    prefix=f"{settings.API_V1_PREFIX}/users",
    tags=["Users"]
)

app.include_router(
    projects.router,
    prefix=f"{settings.API_V1_PREFIX}/projects",
    tags=["Projects"]
)

app.include_router(
    samples.router,
    prefix=f"{settings.API_V1_PREFIX}/samples",
    tags=["Samples"]
)

app.include_router(
    files.router,
    prefix=f"{settings.API_V1_PREFIX}/files",
    tags=["Files"]
)

app.include_router(
    tasks.router,
    prefix=f"{settings.API_V1_PREFIX}/tasks",
    tags=["Tasks"]
)

app.include_router(
    websocket.router,
    prefix=f"{settings.API_V1_PREFIX}",
    tags=["WebSocket"]
)

app.include_router(
    pipelines.router,
    prefix=f"{settings.API_V1_PREFIX}/pipelines",
    tags=["Pipelines"]
)

app.include_router(
    results.router,
    prefix=f"{settings.API_V1_PREFIX}/results",
    tags=["Results"]
)

app.include_router(
    stats.router,
    prefix=f"{settings.API_V1_PREFIX}/stats",
    tags=["Statistics"]
)

app.include_router(
    notifications.router,
    prefix=f"{settings.API_V1_PREFIX}/notifications",
    tags=["Notifications"]
)

app.include_router(
    analytics.router,
    prefix=f"{settings.API_V1_PREFIX}/analytics",
    tags=["Analytics"]
)

app.include_router(
    admin.router,
    prefix=f"{settings.API_V1_PREFIX}/admin",
    tags=["Admin"]
)

app.include_router(
    ai.router,
    prefix=f"{settings.API_V1_PREFIX}/ai",
    tags=["AI Intelligence"]
)

app.include_router(
    data_downloads.router,
    prefix=f"{settings.API_V1_PREFIX}/data-downloads",
    tags=["Data Downloads"]
)

app.include_router(
    vendor_credentials.router,
    prefix=f"{settings.API_V1_PREFIX}/vendor-credentials",
    tags=["Vendor Credentials"]
)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "app.main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.DEBUG
    )
