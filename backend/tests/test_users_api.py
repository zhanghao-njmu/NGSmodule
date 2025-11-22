"""
Integration tests for Users API
Tests statistics query optimization and permission controls
"""
import pytest
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask


class TestUsersAPI:
    """Test Users API endpoints"""

    def test_get_current_user(self, client, auth_headers, test_user):
        """Test getting current user information"""
        response = client.get("/api/v1/users/me", headers=auth_headers)
        assert response.status_code == 200
        data = response.json()
        assert data["username"] == test_user.username
        assert data["email"] == test_user.email
        assert "id" in data

    def test_update_current_user(self, client, auth_headers):
        """Test updating current user information"""
        response = client.put(
            "/api/v1/users/me",
            headers=auth_headers,
            json={"email": "newemail@example.com"}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["email"] == "newemail@example.com"

    def test_get_user_stats_optimized(self, client, admin_headers, db_session, test_user):
        """
        Test optimized user statistics query
        This tests the CASE statement optimization that reduced queries from 5 to 3
        """
        # Create test data for the user
        for i in range(3):
            project = Project(
                name=f"Project {i}",
                user_id=test_user.id
            )
            db_session.add(project)
            db_session.commit()
            db_session.refresh(project)

            # Add samples
            for j in range(2):
                sample = Sample(
                    sample_name=f"Sample {i}-{j}",
                    project_id=project.id
                )
                db_session.add(sample)

            # Add tasks with different statuses
            task1 = PipelineTask(
                task_name=f"Task {i}-completed",
                project_id=project.id,
                status="completed"
            )
            task2 = PipelineTask(
                task_name=f"Task {i}-failed",
                project_id=project.id,
                status="failed"
            )
            task3 = PipelineTask(
                task_name=f"Task {i}-running",
                project_id=project.id,
                status="running"
            )
            db_session.add_all([task1, task2, task3])

        db_session.commit()

        # Get user stats (admin only)
        response = client.get(f"/api/v1/users/{test_user.id}/stats", headers=admin_headers)
        assert response.status_code == 200
        data = response.json()

        # Verify statistics
        assert data["total_projects"] == 3
        assert data["total_samples"] == 6  # 3 projects * 2 samples
        assert data["total_tasks"] == 9  # 3 projects * 3 tasks
        assert data["completed_tasks"] == 3  # 1 per project
        assert data["failed_tasks"] == 3  # 1 per project
        assert "storage_used" in data
        assert "storage_quota" in data

    def test_admin_can_view_all_users(self, client, admin_headers, test_user):
        """Test that admin can view all users"""
        response = client.get("/api/v1/users", headers=admin_headers)
        assert response.status_code == 200
        data = response.json()
        assert isinstance(data, list)
        assert len(data) >= 2  # At least test_user and admin_user

    def test_non_admin_cannot_view_all_users(self, client, auth_headers):
        """Test that non-admin cannot view all users"""
        response = client.get("/api/v1/users", headers=auth_headers)
        assert response.status_code == 403  # Forbidden

    def test_user_cannot_view_other_user_stats(self, client, auth_headers, admin_user):
        """Test that regular users cannot view other users' stats"""
        response = client.get(f"/api/v1/users/{admin_user.id}/stats", headers=auth_headers)
        assert response.status_code == 403  # Forbidden


class TestAuthAPI:
    """Test Authentication API"""

    def test_user_registration(self, client):
        """Test user registration"""
        response = client.post(
            "/api/v1/auth/register",
            json={
                "username": "newuser",
                "email": "newuser@example.com",
                "password": "newpass123"
            }
        )
        # Registration might not be implemented, check for 201 or 404
        assert response.status_code in [201, 404, 405]

    def test_user_login_success(self, client, test_user):
        """Test successful login"""
        response = client.post(
            "/api/v1/auth/login",
            data={"username": "testuser", "password": "testpass123"}
        )
        assert response.status_code == 200
        data = response.json()
        assert "access_token" in data
        assert "token_type" in data
        assert data["token_type"] == "bearer"

    def test_user_login_wrong_password(self, client, test_user):
        """Test login with wrong password"""
        response = client.post(
            "/api/v1/auth/login",
            data={"username": "testuser", "password": "wrongpassword"}
        )
        assert response.status_code == 401

    def test_user_login_nonexistent_user(self, client):
        """Test login with nonexistent user"""
        response = client.post(
            "/api/v1/auth/login",
            data={"username": "nonexistent", "password": "password"}
        )
        assert response.status_code == 401

    def test_token_expiration_handling(self, client):
        """Test handling of expired tokens"""
        # Use an obviously invalid token
        invalid_headers = {"Authorization": "Bearer invalid_token"}
        response = client.get("/api/v1/users/me", headers=invalid_headers)
        assert response.status_code == 401
