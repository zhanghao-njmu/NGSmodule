#!/bin/bash
# NGSmodule Database Restore Script

set -e

# Check if backup file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <backup_file.sql.gz>"
    echo ""
    echo "Available backups:"
    ls -lh "${BACKUP_PATH:-/backups}"/ngsmodule_backup_*.sql.gz 2>/dev/null || echo "No backups found"
    exit 1
fi

BACKUP_FILE="$1"

# Database connection
DB_HOST="${POSTGRES_HOST:-postgres}"
DB_PORT="${POSTGRES_PORT:-5432}"
DB_NAME="${POSTGRES_DB:-ngsmodule}"
DB_USER="${POSTGRES_USER:-ngsmodule}"

echo "==================================="
echo "NGSmodule Database Restore"
echo "==================================="
echo "Timestamp: $(date)"
echo "Database: ${DB_NAME}"
echo "Backup file: ${BACKUP_FILE}"
echo "==================================="

# Confirm restoration
read -p "This will OVERWRITE the current database. Are you sure? (yes/no): " -r
if [[ ! $REPLY =~ ^yes$ ]]; then
    echo "Restore cancelled."
    exit 0
fi

# Check if backup file exists
if [ ! -f "${BACKUP_FILE}" ]; then
    echo "Error: Backup file not found: ${BACKUP_FILE}"
    exit 1
fi

# Decompress if needed
TEMP_FILE="/tmp/restore_temp.sql"
if [[ "${BACKUP_FILE}" == *.gz ]]; then
    echo "Decompressing backup..."
    gunzip -c "${BACKUP_FILE}" > "${TEMP_FILE}"
else
    cp "${BACKUP_FILE}" "${TEMP_FILE}"
fi

# Drop and recreate database
echo "Dropping existing database..."
PGPASSWORD="${POSTGRES_PASSWORD}" psql \
    -h "${DB_HOST}" \
    -p "${DB_PORT}" \
    -U "${DB_USER}" \
    -d postgres \
    -c "DROP DATABASE IF EXISTS ${DB_NAME};"

echo "Creating new database..."
PGPASSWORD="${POSTGRES_PASSWORD}" psql \
    -h "${DB_HOST}" \
    -p "${DB_PORT}" \
    -U "${DB_USER}" \
    -d postgres \
    -c "CREATE DATABASE ${DB_NAME};"

# Restore backup
echo "Restoring backup..."
PGPASSWORD="${POSTGRES_PASSWORD}" pg_restore \
    -h "${DB_HOST}" \
    -p "${DB_PORT}" \
    -U "${DB_USER}" \
    -d "${DB_NAME}" \
    --verbose \
    "${TEMP_FILE}"

# Clean up
rm -f "${TEMP_FILE}"

echo "==================================="
echo "Restore completed successfully"
echo "==================================="
