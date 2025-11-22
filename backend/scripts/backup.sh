#!/bin/bash
# NGSmodule Database Backup Script

set -e

# Configuration
BACKUP_DIR="${BACKUP_PATH:-/backups}"
RETENTION_DAYS="${BACKUP_RETENTION_DAYS:-30}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BACKUP_FILE="${BACKUP_DIR}/ngsmodule_backup_${TIMESTAMP}.sql.gz"

# Database connection
DB_HOST="${POSTGRES_HOST:-postgres}"
DB_PORT="${POSTGRES_PORT:-5432}"
DB_NAME="${POSTGRES_DB:-ngsmodule}"
DB_USER="${POSTGRES_USER:-ngsmodule}"

echo "==================================="
echo "NGSmodule Database Backup"
echo "==================================="
echo "Timestamp: $(date)"
echo "Database: ${DB_NAME}"
echo "Backup file: ${BACKUP_FILE}"
echo "==================================="

# Create backup directory if it doesn't exist
mkdir -p "${BACKUP_DIR}"

# Perform backup
echo "Starting backup..."
PGPASSWORD="${POSTGRES_PASSWORD}" pg_dump \
    -h "${DB_HOST}" \
    -p "${DB_PORT}" \
    -U "${DB_USER}" \
    -d "${DB_NAME}" \
    --verbose \
    --format=custom \
    --file="${BACKUP_FILE%.gz}"

# Compress backup
echo "Compressing backup..."
gzip "${BACKUP_FILE%.gz}"

# Calculate backup size
BACKUP_SIZE=$(du -h "${BACKUP_FILE}" | cut -f1)
echo "Backup completed: ${BACKUP_SIZE}"

# Remove old backups
echo "Cleaning up old backups (older than ${RETENTION_DAYS} days)..."
find "${BACKUP_DIR}" -name "ngsmodule_backup_*.sql.gz" -type f -mtime +${RETENTION_DAYS} -delete

echo "==================================="
echo "Backup completed successfully"
echo "==================================="

# List recent backups
echo ""
echo "Recent backups:"
ls -lh "${BACKUP_DIR}"/ngsmodule_backup_*.sql.gz | tail -5
