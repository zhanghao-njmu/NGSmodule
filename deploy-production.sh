#!/bin/bash
# ================================================================================
# NGSmodule Production Deployment Script
# ================================================================================
# This script automates the production deployment process
# Usage: ./deploy-production.sh [command]
#
# Commands:
#   setup     - Initial setup (create directories, generate secrets)
#   build     - Build Docker images
#   start     - Start all services
#   stop      - Stop all services
#   restart   - Restart all services
#   logs      - View logs
#   status    - Check service status
#   backup    - Backup database
#   restore   - Restore database from backup
#   update    - Pull latest code and rebuild
#   clean     - Clean up (remove containers, volumes, images)
# ================================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PROJECT_NAME="ngsmodule"
COMPOSE_FILE="docker-compose.prod.yml"
ENV_FILE=".env"
ENV_TEMPLATE=".env.production"
BACKUP_DIR="./backups"

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_requirements() {
    log_info "Checking requirements..."

    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed. Please install Docker first."
        exit 1
    fi

    if ! command -v docker-compose &> /dev/null; then
        log_error "Docker Compose is not installed. Please install Docker Compose first."
        exit 1
    fi

    log_success "All requirements met"
}

setup() {
    log_info "Starting initial setup..."

    # Check if .env already exists
    if [ -f "$ENV_FILE" ]; then
        log_warning ".env file already exists. Skipping creation."
        log_info "If you want to regenerate, delete .env and run setup again."
    else
        # Copy template
        if [ -f "$ENV_TEMPLATE" ]; then
            cp "$ENV_TEMPLATE" "$ENV_FILE"
            log_success "Created .env from template"
        else
            log_error ".env.production template not found"
            exit 1
        fi

        # Generate secrets
        log_info "Generating secure random secrets..."
        SECRET_KEY=$(openssl rand -hex 32)
        JWT_SECRET_KEY=$(openssl rand -hex 32)
        POSTGRES_PASSWORD=$(openssl rand -base64 32 | tr -d "=+/" | cut -c1-32)
        REDIS_PASSWORD=$(openssl rand -base64 32 | tr -d "=+/" | cut -c1-32)
        MINIO_PASSWORD=$(openssl rand -base64 32 | tr -d "=+/" | cut -c1-32)
        FLOWER_PASSWORD=$(openssl rand -base64 16 | tr -d "=+/" | cut -c1-16)

        # Update .env file
        sed -i "s/CHANGE_THIS_TO_64_CHAR_HEX_STRING_use_openssl_rand_hex_32/${SECRET_KEY}/g" "$ENV_FILE"
        sed -i "0,/CHANGE_THIS_TO_64_CHAR_HEX_STRING_use_openssl_rand_hex_32/s//${JWT_SECRET_KEY}/" "$ENV_FILE"
        sed -i "s/CHANGE_THIS_TO_STRONG_PASSWORD_min_32_chars/${POSTGRES_PASSWORD}/g" "$ENV_FILE"
        sed -i "s/CHANGE_THIS_TO_STRONG_REDIS_PASSWORD_min_32_chars/${REDIS_PASSWORD}/g" "$ENV_FILE"
        sed -i "s/CHANGE_THIS_TO_STRONG_MINIO_PASSWORD_min_32_chars/${MINIO_PASSWORD}/g" "$ENV_FILE"
        sed -i "s/CHANGE_THIS_FLOWER_PASSWORD_min_16_chars/${FLOWER_PASSWORD}/g" "$ENV_FILE"

        log_success "Generated and saved secure secrets"
        log_warning "IMPORTANT: Review and update .env file with your production values:"
        log_warning "  - Update domain names (yourdomain.com)"
        log_warning "  - Configure email settings (if needed)"
        log_warning "  - Set MINIO_ROOT_USER"
        log_warning "  - Review other settings"
    fi

    # Create necessary directories
    log_info "Creating required directories..."
    mkdir -p backups/postgres
    mkdir -p ssl
    mkdir -p logs

    # Set permissions
    chmod 600 "$ENV_FILE"
    chmod 700 backups ssl logs

    log_success "Directory structure created"

    # SSL certificate setup
    if [ ! -f "ssl/cert.pem" ] || [ ! -f "ssl/key.pem" ]; then
        log_warning "SSL certificates not found in ./ssl/"
        log_info "You can:"
        log_info "  1. Use Let's Encrypt: certbot certonly --standalone -d yourdomain.com"
        log_info "  2. Generate self-signed cert: openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout ssl/key.pem -out ssl/cert.pem"
        log_info "  3. Copy existing certificates to ./ssl/cert.pem and ./ssl/key.pem"
    fi

    log_success "Setup completed!"
    log_info "Next steps:"
    log_info "  1. Edit .env and update production values"
    log_info "  2. Set up SSL certificates in ./ssl/"
    log_info "  3. Run: ./deploy-production.sh build"
}

build() {
    log_info "Building Docker images..."

    # Pull base images
    docker-compose -f "$COMPOSE_FILE" pull

    # Build application images
    docker-compose -f "$COMPOSE_FILE" build --no-cache

    log_success "Docker images built successfully"
}

start() {
    log_info "Starting services..."

    docker-compose -f "$COMPOSE_FILE" up -d

    log_success "Services started"
    log_info "Waiting for services to be healthy..."
    sleep 10

    status
}

stop() {
    log_info "Stopping services..."

    docker-compose -f "$COMPOSE_FILE" down

    log_success "Services stopped"
}

restart() {
    log_info "Restarting services..."

    stop
    sleep 5
    start
}

logs() {
    log_info "Showing logs (Ctrl+C to exit)..."

    docker-compose -f "$COMPOSE_FILE" logs -f --tail=100
}

status() {
    log_info "Service status:"
    echo ""

    docker-compose -f "$COMPOSE_FILE" ps

    echo ""
    log_info "Health checks:"

    # Check backend health
    if curl -sf http://localhost:8000/health > /dev/null 2>&1; then
        log_success "Backend API: Healthy"
    else
        log_error "Backend API: Unhealthy or not responding"
    fi

    # Check frontend health
    if curl -sf http://localhost/health > /dev/null 2>&1; then
        log_success "Frontend: Healthy"
    else
        log_error "Frontend: Unhealthy or not responding"
    fi

    # Check nginx status
    if curl -sf http://localhost:8080/nginx-status > /dev/null 2>&1; then
        log_success "Nginx: Healthy"
    else
        log_warning "Nginx: Status endpoint not accessible"
    fi
}

backup() {
    log_info "Backing up PostgreSQL database..."

    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    BACKUP_FILE="${BACKUP_DIR}/postgres/backup_${TIMESTAMP}.sql"

    docker-compose -f "$COMPOSE_FILE" exec -T postgres pg_dump -U ngsmodule ngsmodule > "$BACKUP_FILE"

    # Compress backup
    gzip "$BACKUP_FILE"

    log_success "Database backed up to: ${BACKUP_FILE}.gz"

    # Clean old backups (keep last 30 days)
    find "${BACKUP_DIR}/postgres" -name "backup_*.sql.gz" -mtime +30 -delete
    log_info "Old backups cleaned up (retention: 30 days)"
}

restore() {
    if [ -z "$1" ]; then
        log_error "Please specify backup file to restore"
        log_info "Usage: ./deploy-production.sh restore <backup_file>"
        log_info "Available backups:"
        ls -lh "${BACKUP_DIR}/postgres/"
        exit 1
    fi

    BACKUP_FILE="$1"

    if [ ! -f "$BACKUP_FILE" ]; then
        log_error "Backup file not found: $BACKUP_FILE"
        exit 1
    fi

    log_warning "This will restore database from: $BACKUP_FILE"
    read -p "Are you sure? (yes/no): " confirm

    if [ "$confirm" != "yes" ]; then
        log_info "Restore cancelled"
        exit 0
    fi

    log_info "Restoring database..."

    # Decompress if needed
    if [[ "$BACKUP_FILE" == *.gz ]]; then
        gunzip -c "$BACKUP_FILE" | docker-compose -f "$COMPOSE_FILE" exec -T postgres psql -U ngsmodule ngsmodule
    else
        docker-compose -f "$COMPOSE_FILE" exec -T postgres psql -U ngsmodule ngsmodule < "$BACKUP_FILE"
    fi

    log_success "Database restored successfully"
}

update() {
    log_info "Updating application..."

    # Pull latest code
    log_info "Pulling latest code from git..."
    git pull

    # Rebuild images
    build

    # Restart services
    restart

    log_success "Application updated successfully"
}

clean() {
    log_warning "This will remove all containers, volumes, and images"
    read -p "Are you sure? (yes/no): " confirm

    if [ "$confirm" != "yes" ]; then
        log_info "Clean cancelled"
        exit 0
    fi

    log_info "Cleaning up..."

    # Stop and remove containers
    docker-compose -f "$COMPOSE_FILE" down -v

    # Remove images
    docker-compose -f "$COMPOSE_FILE" down --rmi all

    log_success "Cleanup completed"
}

# Main script
case "${1:-}" in
    setup)
        check_requirements
        setup
        ;;
    build)
        check_requirements
        build
        ;;
    start)
        check_requirements
        start
        ;;
    stop)
        check_requirements
        stop
        ;;
    restart)
        check_requirements
        restart
        ;;
    logs)
        check_requirements
        logs
        ;;
    status)
        check_requirements
        status
        ;;
    backup)
        check_requirements
        backup
        ;;
    restore)
        check_requirements
        restore "$2"
        ;;
    update)
        check_requirements
        update
        ;;
    clean)
        check_requirements
        clean
        ;;
    *)
        echo "NGSmodule Production Deployment Script"
        echo ""
        echo "Usage: $0 [command]"
        echo ""
        echo "Commands:"
        echo "  setup     - Initial setup (create directories, generate secrets)"
        echo "  build     - Build Docker images"
        echo "  start     - Start all services"
        echo "  stop      - Stop all services"
        echo "  restart   - Restart all services"
        echo "  logs      - View logs"
        echo "  status    - Check service status"
        echo "  backup    - Backup database"
        echo "  restore   - Restore database from backup"
        echo "  update    - Pull latest code and rebuild"
        echo "  clean     - Clean up (remove containers, volumes, images)"
        echo ""
        echo "Example:"
        echo "  $0 setup      # First time setup"
        echo "  $0 build      # Build images"
        echo "  $0 start      # Start services"
        echo "  $0 status     # Check health"
        echo ""
        exit 1
        ;;
esac
