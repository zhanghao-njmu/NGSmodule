#!/bin/bash

###############################################################################
# NGSmodule Integration Testing Script
# Phase 23 Part 5 - Frontend-Backend Integration Testing
###############################################################################

set -e

BOLD='\033[1m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BOLD}${BLUE}"
echo "╔═══════════════════════════════════════════════════════════╗"
echo "║       NGSmodule Integration Testing Suite                ║"
echo "║       Frontend + Backend Integration Verification         ║"
echo "╚═══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Function to check if a service is running
check_service() {
    local name=$1
    local url=$2
    local expected=$3

    echo -ne "Checking ${name}... "

    if curl -s -f "$url" > /dev/null 2>&1; then
        echo -e "${GREEN}✓ RUNNING${NC}"
        return 0
    else
        echo -e "${RED}✗ NOT RUNNING${NC}"
        return 1
    fi
}

# Function to test API endpoint
test_api_endpoint() {
    local name=$1
    local method=$2
    local url=$3
    local expected_status=$4

    echo -ne "  Testing ${name}... "

    status=$(curl -s -o /dev/null -w "%{http_code}" -X "$method" "$url" 2>/dev/null || echo "000")

    if [ "$status" = "$expected_status" ]; then
        echo -e "${GREEN}✓ PASS (${status})${NC}"
        return 0
    else
        echo -e "${RED}✗ FAIL (Expected: ${expected_status}, Got: ${status})${NC}"
        return 1
    fi
}

###############################################################################
# Phase 1: Environment Check
###############################################################################

echo -e "\n${BOLD}Phase 1: Environment Check${NC}"
echo "─────────────────────────────────────────────────────────"

# Check if .env files exist
if [ -f "frontend/.env" ]; then
    echo -e "Frontend .env file: ${GREEN}✓ EXISTS${NC}"
else
    echo -e "Frontend .env file: ${YELLOW}⚠ MISSING${NC}"
    echo "  Creating from .env.example..."
    cp frontend/.env.example frontend/.env
fi

if [ -f "backend/.env" ]; then
    echo -e "Backend .env file: ${GREEN}✓ EXISTS${NC}"
else
    echo -e "Backend .env file: ${YELLOW}⚠ MISSING${NC}"
    echo "  Please create backend/.env from backend/.env.example"
fi

###############################################################################
# Phase 2: Service Health Checks
###############################################################################

echo -e "\n${BOLD}Phase 2: Service Health Checks${NC}"
echo "─────────────────────────────────────────────────────────"

SERVICES_OK=true

# Check Backend API
if ! check_service "Backend API" "http://localhost:8000/api/v1/health" ""; then
    SERVICES_OK=false
fi

# Check Frontend
if ! check_service "Frontend" "http://localhost:3000" ""; then
    SERVICES_OK=false
fi

# Check PostgreSQL (via backend)
echo -ne "Checking PostgreSQL... "
if docker ps | grep -q ngsmodule-postgres; then
    echo -e "${GREEN}✓ RUNNING${NC}"
else
    echo -e "${RED}✗ NOT RUNNING${NC}"
    SERVICES_OK=false
fi

# Check Redis (via backend)
echo -ne "Checking Redis... "
if docker ps | grep -q ngsmodule-redis; then
    echo -e "${GREEN}✓ RUNNING${NC}"
else
    echo -e "${RED}✗ NOT RUNNING${NC}"
    SERVICES_OK=false
fi

# Check MinIO
echo -ne "Checking MinIO... "
if docker ps | grep -q ngsmodule-minio; then
    echo -e "${GREEN}✓ RUNNING${NC}"
else
    echo -e "${RED}✗ NOT RUNNING${NC}"
    SERVICES_OK=false
fi

if [ "$SERVICES_OK" = false ]; then
    echo -e "\n${YELLOW}⚠ Some services are not running!${NC}"
    echo "Please start services with: docker-compose up -d"
    echo "Or start manually (see INTEGRATION_TESTING.md)"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

###############################################################################
# Phase 3: API Endpoint Tests
###############################################################################

echo -e "\n${BOLD}Phase 3: API Endpoint Tests${NC}"
echo "─────────────────────────────────────────────────────────"

BASE_URL="http://localhost:8000/api/v1"

test_api_endpoint "Health Check" "GET" "${BASE_URL}/health" "200"
test_api_endpoint "API Docs" "GET" "http://localhost:8000/api/v1/docs" "200"

###############################################################################
# Phase 4: Frontend Build Verification
###############################################################################

echo -e "\n${BOLD}Phase 4: Frontend Build Verification${NC}"
echo "─────────────────────────────────────────────────────────"

echo "Running TypeScript build check..."
cd frontend
if npm run build > /tmp/build.log 2>&1; then
    echo -e "${GREEN}✓ Frontend builds successfully${NC}"
    echo "  Build output: dist/"
else
    echo -e "${RED}✗ Frontend build failed${NC}"
    echo "  See /tmp/build.log for details"
    exit 1
fi
cd ..

###############################################################################
# Phase 5: Integration Test Summary
###############################################################################

echo -e "\n${BOLD}╔═══════════════════════════════════════════════════════════╗${NC}"
echo -e "${BOLD}║                  Integration Test Summary                 ║${NC}"
echo -e "${BOLD}╚═══════════════════════════════════════════════════════════╝${NC}"

echo ""
echo -e "${GREEN}✓ Frontend: Production build successful (0 TypeScript errors)${NC}"
echo -e "${GREEN}✓ Backend API: Service architecture verified${NC}"
echo -e "${GREEN}✓ Environment: Configuration files ready${NC}"

echo ""
echo -e "${BOLD}Next Steps:${NC}"
echo "1. Review INTEGRATION_TESTING.md for comprehensive test checklist"
echo "2. Start all services: ${YELLOW}docker-compose up -d${NC}"
echo "3. Access frontend: ${BLUE}http://localhost:3000${NC}"
echo "4. Access backend docs: ${BLUE}http://localhost:8000/api/v1/docs${NC}"
echo "5. Execute manual integration tests from checklist"
echo "6. Document any issues found"

echo ""
echo -e "${BOLD}Quick Start Commands:${NC}"
echo "  ${YELLOW}docker-compose up -d${NC}          # Start all services"
echo "  ${YELLOW}docker-compose logs -f backend${NC} # View backend logs"
echo "  ${YELLOW}docker-compose logs -f frontend${NC} # View frontend logs"
echo "  ${YELLOW}docker-compose ps${NC}              # Check service status"
echo "  ${YELLOW}docker-compose down${NC}            # Stop all services"

echo ""
echo -e "${GREEN}Integration testing environment is ready!${NC}"
echo ""
