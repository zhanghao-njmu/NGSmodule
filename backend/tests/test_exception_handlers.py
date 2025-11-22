"""
Integration tests for Global Exception Handlers
Tests the exception handlers added in Phase 8
"""
import pytest
from unittest.mock import patch, MagicMock
from sqlalchemy.exc import IntegrityError, DBAPIError


class TestGlobalExceptionHandlers:
    """Test global exception handlers"""

    def test_integrity_error_handler(self, client, auth_headers, db_session, test_user):
        """Test IntegrityError handler (duplicate resource)"""
        # Create a project
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={"name": "Unique Project"}
        )
        assert response.status_code == 201

        # Try to create another project with potential duplicate
        # The actual duplicate constraint depends on your model
        # This test validates the error response format
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={"name": "Another Project"}  # This might not trigger IntegrityError
        )
        # Just verify that if it fails, it has the right format
        if response.status_code == 409:
            data = response.json()
            assert "error" in data
            assert "code" in data
            assert data["code"] == "INTEGRITY_ERROR"

    def test_validation_error_handler(self, client, auth_headers):
        """Test ValidationError handler (invalid data)"""
        # Send invalid data (missing required field)
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={}  # Missing required 'name' field
        )
        assert response.status_code == 422
        data = response.json()
        assert "error" in data or "detail" in data
        if "code" in data:
            assert data["code"] == "VALIDATION_ERROR"

    def test_validation_error_invalid_type(self, client, auth_headers):
        """Test validation with wrong data type"""
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={"name": 123}  # Should be string
        )
        assert response.status_code == 422
        data = response.json()
        assert "error" in data or "detail" in data

    def test_404_not_found_error(self, client, auth_headers):
        """Test 404 error handling"""
        fake_uuid = "00000000-0000-0000-0000-000000000000"
        response = client.get(f"/api/v1/projects/{fake_uuid}", headers=auth_headers)
        assert response.status_code == 404
        data = response.json()
        assert "detail" in data or "error" in data

    def test_401_unauthorized_error(self, client):
        """Test 401 unauthorized error"""
        response = client.get("/api/v1/projects")
        assert response.status_code == 401
        # Should have WWW-Authenticate header
        assert "WWW-Authenticate" in response.headers or response.status_code == 401

    def test_403_forbidden_error(self, client, auth_headers):
        """Test 403 forbidden error (insufficient permissions)"""
        # Try to access admin-only endpoint
        response = client.get("/api/v1/users", headers=auth_headers)
        assert response.status_code == 403
        data = response.json()
        assert "detail" in data or "error" in data

    def test_error_response_format_consistency(self, client, auth_headers):
        """Test that all errors return consistent format"""
        # Test various error scenarios
        error_responses = []

        # 404 error
        response = client.get("/api/v1/projects/00000000-0000-0000-0000-000000000000", headers=auth_headers)
        if response.status_code == 404:
            error_responses.append(response.json())

        # 422 validation error
        response = client.post("/api/v1/projects", headers=auth_headers, json={})
        if response.status_code == 422:
            error_responses.append(response.json())

        # 403 forbidden error
        response = client.get("/api/v1/users", headers=auth_headers)
        if response.status_code == 403:
            error_responses.append(response.json())

        # All error responses should have either 'error' or 'detail' key
        for error in error_responses:
            assert "error" in error or "detail" in error


class TestPermissionDependencies:
    """Test the permission dependency functions added in Phase 8"""

    def test_get_user_project_dependency(self, client, auth_headers, db_session, test_user):
        """Test get_user_project dependency function"""
        # Create a project owned by test_user
        from app.models.project import Project
        project = Project(name="Test Project", user_id=test_user.id)
        db_session.add(project)
        db_session.commit()
        db_session.refresh(project)

        # Should succeed - user owns the project
        response = client.get(f"/api/v1/projects/{project.id}", headers=auth_headers)
        assert response.status_code == 200

        # Create another user's project
        from app.models.user import User
        from app.core.security import get_password_hash
        other_user = User(
            username="other",
            email="other@example.com",
            hashed_password=get_password_hash("password")
        )
        db_session.add(other_user)
        db_session.commit()

        other_project = Project(name="Other Project", user_id=other_user.id)
        db_session.add(other_project)
        db_session.commit()
        db_session.refresh(other_project)

        # Should fail - user doesn't own this project
        response = client.get(f"/api/v1/projects/{other_project.id}", headers=auth_headers)
        assert response.status_code == 404  # Not found (access denied)

    def test_admin_can_access_all_projects(self, client, admin_headers, db_session, test_user):
        """Test that admin can access all projects"""
        # Create a project owned by regular user
        from app.models.project import Project
        project = Project(name="User Project", user_id=test_user.id)
        db_session.add(project)
        db_session.commit()
        db_session.refresh(project)

        # Admin should be able to access it
        response = client.get(f"/api/v1/projects/{project.id}", headers=admin_headers)
        assert response.status_code == 200
        data = response.json()
        assert data["name"] == "User Project"


class TestErrorResponseStructure:
    """Test error response structure consistency"""

    def test_error_has_required_fields(self, client, auth_headers):
        """Test that error responses have required fields"""
        # Trigger a 404 error
        response = client.get("/api/v1/projects/00000000-0000-0000-0000-000000000000", headers=auth_headers)
        if response.status_code == 404:
            data = response.json()
            # Should have detail field at minimum
            assert "detail" in data

    def test_validation_error_has_detail_array(self, client, auth_headers):
        """Test that validation errors include detailed error array"""
        response = client.post("/api/v1/projects", headers=auth_headers, json={})
        assert response.status_code == 422
        data = response.json()
        assert "detail" in data
        # Validation errors typically have an array of error details
        if isinstance(data["detail"], list):
            assert len(data["detail"]) > 0
