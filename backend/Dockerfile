# Multi-stage build for production

# Stage 1: Builder
FROM python:3.11-slim as builder

WORKDIR /app

# Install build dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Build Python deps inside a virtualenv at /opt/venv. This puts both
# site-packages and the bin/ scripts (uvicorn, celery, alembic, ...) at
# a path that the non-root runtime user can read in stage 2 — installing
# with `pip install --user` would put them under /root/.local, which is
# 0700 root-owned and crashes the container with EACCES on uvicorn at boot.
RUN python -m venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH

# Copy requirements and install
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Stage 2: Runtime
FROM python:3.11-slim

WORKDIR /app

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libpq5 \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy the prepared venv from the builder. /opt/venv is world-readable,
# so the non-root user we add below can invoke uvicorn/celery/alembic.
COPY --from=builder /opt/venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /app/storage /app/logs

# Create non-root user
RUN useradd -m -u 1000 ngsmodule && \
    chown -R ngsmodule:ngsmodule /app

# Switch to non-root user
USER ngsmodule

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Run application
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "4"]
