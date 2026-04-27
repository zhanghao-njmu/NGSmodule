#!/bin/bash

###############################################################################
# NGSmodule 前后端 API 集成测试脚本
#
# 功能：自动测试所有后端 API 端点，验证前端集成
# 用法：./test-api-integration.sh [BASE_URL]
# 例如：./test-api-integration.sh http://localhost:8000
###############################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration
API_URL="${1:-http://localhost:8000}"
API_PREFIX="/api/v1"
TEST_USERNAME="integration_test_$(date +%s)"
TEST_PASSWORD="TestPass123!"
TEST_EMAIL="${TEST_USERNAME}@test.com"

# Counters
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0
SKIPPED_TESTS=0

# Variables to store data between tests
ACCESS_TOKEN=""
USER_ID=""
PROJECT_ID=""
SAMPLE_ID=""

###############################################################################
# Helper Functions
###############################################################################

print_header() {
    echo ""
    echo -e "${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║ $1${NC}"
    echo -e "${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
}

print_section() {
    echo ""
    echo -e "${BLUE}── $1 ──${NC}"
}

print_success() {
    echo -e "  ${GREEN}✓${NC} $1"
    PASSED_TESTS=$((PASSED_TESTS + 1))
}

print_failure() {
    echo -e "  ${RED}✗${NC} $1"
    if [ -n "$2" ]; then
        echo -e "    ${RED}Error: $2${NC}"
    fi
    FAILED_TESTS=$((FAILED_TESTS + 1))
}

print_skip() {
    echo -e "  ${YELLOW}○${NC} $1 (skipped)"
    SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
}

print_info() {
    echo -e "  ${YELLOW}ℹ${NC} $1"
}

# Test endpoint - returns 0 if success, 1 if failure
test_endpoint() {
    local method=$1
    local endpoint=$2
    local expected_status=$3
    local description=$4
    local data=$5
    local auth=$6

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    local curl_args=("-s" "-o" "/tmp/response.json" "-w" "%{http_code}" "-X" "$method")

    if [ "$auth" = "true" ] && [ -n "$ACCESS_TOKEN" ]; then
        curl_args+=("-H" "Authorization: Bearer $ACCESS_TOKEN")
    fi

    if [ -n "$data" ]; then
        curl_args+=("-H" "Content-Type: application/json" "-d" "$data")
    fi

    local url="${API_URL}${API_PREFIX}${endpoint}"
    local status=$(curl "${curl_args[@]}" "$url" 2>/dev/null)

    if [ "$status" = "$expected_status" ]; then
        print_success "$description (HTTP $status)"
        return 0
    else
        local error_msg=""
        if [ -f /tmp/response.json ]; then
            error_msg=$(cat /tmp/response.json 2>/dev/null | head -c 200)
        fi
        print_failure "$description (Expected $expected_status, got $status)" "$error_msg"
        return 1
    fi
}

# Extract value from JSON response
extract_json() {
    local key=$1
    if [ -f /tmp/response.json ]; then
        grep -o "\"$key\":\"[^\"]*\"" /tmp/response.json | head -1 | sed "s/\"$key\":\"\([^\"]*\)\"/\1/"
    fi
}

###############################################################################
# Test Suites
###############################################################################

test_health_check() {
    print_section "1. Health Check"

    test_endpoint "GET" "/../health" "200" "Server health check"
    test_endpoint "GET" "/openapi.json" "200" "OpenAPI specification accessible"
}

test_authentication() {
    print_section "2. Authentication Flow"

    # Register new user
    local register_data="{\"username\":\"$TEST_USERNAME\",\"email\":\"$TEST_EMAIL\",\"password\":\"$TEST_PASSWORD\",\"full_name\":\"Integration Test User\"}"
    if test_endpoint "POST" "/auth/register" "201" "User registration" "$register_data"; then
        USER_ID=$(extract_json "id")
        print_info "Created user with ID: $USER_ID"
    fi

    # Login
    local login_data="{\"username\":\"$TEST_USERNAME\",\"password\":\"$TEST_PASSWORD\"}"
    if test_endpoint "POST" "/auth/login" "200" "User login" "$login_data"; then
        ACCESS_TOKEN=$(grep -o '"access_token":"[^"]*"' /tmp/response.json | head -1 | sed 's/"access_token":"\([^"]*\)"/\1/')
        if [ -n "$ACCESS_TOKEN" ]; then
            print_info "Got access token (length: ${#ACCESS_TOKEN})"
        else
            print_failure "Could not extract access token from response"
        fi
    fi

    # Get current user info
    if [ -n "$ACCESS_TOKEN" ]; then
        test_endpoint "GET" "/auth/me" "200" "Get current user info" "" "true"
    else
        print_skip "Get current user info"
    fi
}

test_users() {
    print_section "3. Users API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All user tests (no token)"
        return
    fi

    test_endpoint "GET" "/users/me" "200" "Get my profile" "" "true"
    test_endpoint "GET" "/users/me/storage" "200" "Get my storage info" "" "true"
}

test_projects() {
    print_section "4. Projects API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All project tests (no token)"
        return
    fi

    # List projects
    test_endpoint "GET" "/projects" "200" "List projects" "" "true"

    # Create project
    local project_data='{"name":"Integration Test Project","description":"Created by integration test"}'
    if test_endpoint "POST" "/projects" "201" "Create project" "$project_data" "true"; then
        PROJECT_ID=$(extract_json "id")
        if [ -n "$PROJECT_ID" ]; then
            print_info "Created project with ID: $PROJECT_ID"

            # Get project
            test_endpoint "GET" "/projects/$PROJECT_ID" "200" "Get project by ID" "" "true"

            # Update project
            local update_data='{"description":"Updated description"}'
            test_endpoint "PUT" "/projects/$PROJECT_ID" "200" "Update project" "$update_data" "true"
        fi
    fi
}

test_samples() {
    print_section "5. Samples API"

    if [ -z "$ACCESS_TOKEN" ] || [ -z "$PROJECT_ID" ]; then
        print_skip "All sample tests (no token or project)"
        return
    fi

    # List samples
    test_endpoint "GET" "/samples?project_id=$PROJECT_ID" "200" "List project samples" "" "true"
}

test_files() {
    print_section "6. Files API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All file tests (no token)"
        return
    fi

    test_endpoint "GET" "/files" "200" "List files" "" "true"
}

test_tasks() {
    print_section "7. Tasks API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All task tests (no token)"
        return
    fi

    test_endpoint "GET" "/tasks" "200" "List tasks" "" "true"
}

test_pipelines() {
    print_section "8. Pipelines API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All pipeline tests (no token)"
        return
    fi

    test_endpoint "GET" "/pipelines" "200" "List pipelines" "" "true"
}

test_results() {
    print_section "9. Results API"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All result tests (no token)"
        return
    fi

    test_endpoint "GET" "/results" "200" "List results" "" "true"
}

test_stats() {
    print_section "10. Stats API (NEW)"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All stats tests (no token)"
        return
    fi

    test_endpoint "GET" "/stats/summary" "200" "Get stats summary" "" "true"
    test_endpoint "GET" "/stats/projects" "200" "Get project stats" "" "true"
    test_endpoint "GET" "/stats/samples" "200" "Get sample stats" "" "true"
    test_endpoint "GET" "/stats/tasks" "200" "Get task stats" "" "true"
    test_endpoint "GET" "/stats/files" "200" "Get file stats" "" "true"
    test_endpoint "GET" "/stats/storage" "200" "Get storage stats" "" "true"
    test_endpoint "GET" "/stats/pipelines" "200" "Get pipeline stats" "" "true"
    test_endpoint "GET" "/stats/quick" "200" "Get quick stats (dashboard)" "" "true"
    test_endpoint "GET" "/stats/trends/tasks?period=daily&days=7" "200" "Get task trends" "" "true"
}

test_notifications() {
    print_section "11. Notifications API (NEW)"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All notification tests (no token)"
        return
    fi

    test_endpoint "GET" "/notifications" "200" "Get notifications list" "" "true"
    test_endpoint "GET" "/notifications/unread/count" "200" "Get unread count" "" "true"
    test_endpoint "GET" "/notifications/settings/current" "200" "Get notification settings" "" "true"

    # Update settings
    local settings_data='{"email_enabled":true,"app_enabled":true,"push_enabled":false}'
    test_endpoint "PUT" "/notifications/settings/current" "200" "Update notification settings" "$settings_data" "true"
}

test_analytics() {
    print_section "12. Analytics API (NEW)"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All analytics tests (no token)"
        return
    fi

    test_endpoint "GET" "/analytics/dashboard" "200" "Get dashboard analytics" "" "true"
    test_endpoint "GET" "/analytics/timeseries/tasks?time_range=week" "200" "Get tasks time series" "" "true"
    test_endpoint "GET" "/analytics/trends/tasks?time_range=month" "200" "Get task trends" "" "true"
    test_endpoint "GET" "/analytics/storage" "200" "Get storage analytics" "" "true"
    test_endpoint "GET" "/analytics/pipelines/performance" "200" "Get pipeline performance" "" "true"
    test_endpoint "GET" "/analytics/tasks/execution-trend?time_range=week" "200" "Get task execution trend" "" "true"
    test_endpoint "GET" "/analytics/health" "200" "Analytics health check"
}

test_admin() {
    print_section "13. Admin API (NEW - Admin only)"

    if [ -z "$ACCESS_TOKEN" ]; then
        print_skip "All admin tests (no token)"
        return
    fi

    # These will return 403 for non-admin users (expected)
    print_info "Note: These endpoints require admin role. 403 responses are expected for regular users."

    # Test that endpoints exist (403 is OK for non-admin)
    local status=$(curl -s -o /tmp/response.json -w "%{http_code}" -X GET \
        -H "Authorization: Bearer $ACCESS_TOKEN" \
        "${API_URL}${API_PREFIX}/admin/users" 2>/dev/null)

    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    if [ "$status" = "200" ] || [ "$status" = "403" ]; then
        print_success "Admin /users endpoint exists (HTTP $status)"
    else
        print_failure "Admin /users endpoint test (got $status)" ""
    fi

    # Check other admin endpoints
    for endpoint in "/admin/config" "/admin/logs" "/admin/system/health" "/admin/system/stats"; do
        local status=$(curl -s -o /tmp/response.json -w "%{http_code}" -X GET \
            -H "Authorization: Bearer $ACCESS_TOKEN" \
            "${API_URL}${API_PREFIX}${endpoint}" 2>/dev/null)

        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        if [ "$status" = "200" ] || [ "$status" = "403" ]; then
            print_success "Admin ${endpoint} exists (HTTP $status)"
        else
            print_failure "Admin ${endpoint} test (got $status)" ""
        fi
    done
}

test_cleanup() {
    print_section "14. Cleanup"

    # Delete test project
    if [ -n "$ACCESS_TOKEN" ] && [ -n "$PROJECT_ID" ]; then
        test_endpoint "DELETE" "/projects/$PROJECT_ID" "204" "Delete test project" "" "true"
    fi

    print_info "Test user '$TEST_USERNAME' was created. Manual cleanup may be needed."
}

###############################################################################
# Main Execution
###############################################################################

main() {
    print_header "NGSmodule API Integration Tests"
    echo "API URL: ${API_URL}${API_PREFIX}"
    echo "Timestamp: $(date)"
    echo ""

    # Check if server is reachable
    if ! curl -s -f -o /dev/null "${API_URL}/health"; then
        echo -e "${RED}❌ ERROR: Cannot reach server at ${API_URL}${NC}"
        echo "Please ensure the backend is running:"
        echo "  cd backend && uvicorn app.main:app --reload"
        exit 1
    fi

    # Run all test suites
    test_health_check
    test_authentication
    test_users
    test_projects
    test_samples
    test_files
    test_tasks
    test_pipelines
    test_results
    test_stats           # NEW Phase 28
    test_notifications   # NEW Phase 28
    test_analytics       # NEW Phase 28
    test_admin           # NEW Phase 28
    test_cleanup

    # Print summary
    print_header "Test Summary"
    echo ""
    echo "  Total Tests:   $TOTAL_TESTS"
    echo -e "  ${GREEN}Passed:        $PASSED_TESTS${NC}"
    echo -e "  ${RED}Failed:        $FAILED_TESTS${NC}"
    echo -e "  ${YELLOW}Skipped:       $SKIPPED_TESTS${NC}"
    echo ""

    # Calculate pass rate
    if [ $TOTAL_TESTS -gt 0 ]; then
        local pass_rate=$(echo "scale=2; $PASSED_TESTS * 100 / $TOTAL_TESTS" | bc 2>/dev/null || echo "N/A")
        echo "  Pass Rate:     ${pass_rate}%"
    fi

    echo ""

    if [ $FAILED_TESTS -eq 0 ]; then
        echo -e "${GREEN}✅ All tests passed!${NC}"
        exit 0
    else
        echo -e "${RED}❌ Some tests failed. Please review the output above.${NC}"
        exit 1
    fi
}

# Run main function
main
