"""
Integration tests for Projects API
Tests N+1 query optimization, permissions, and response format
"""
import pytest
from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask


class TestProjectsAPI:
    """Test Projects API endpoints"""

    def test_create_project(self, client, auth_headers):
        """Test creating a new project"""
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={
                "name": "Test Project",
                "description": "Test Description",
                "project_type": "WGS"
            }
        )
        assert response.status_code == 201
        data = response.json()
        assert data["name"] == "Test Project"
        assert data["description"] == "Test Description"
        assert data["project_type"] == "WGS"
        assert "id" in data
        assert "created_at" in data

    def test_list_projects_empty(self, client, auth_headers):
        """Test listing projects when none exist"""
        response = client.get("/api/v1/projects", headers=auth_headers)
        assert response.status_code == 200
        data = response.json()
        assert data["total"] == 0
        assert data["items"] == []
        assert "page" in data
        assert "page_size" in data

    def test_list_projects_with_counts(self, client, auth_headers, db_session, test_user):
        """
        Test N+1 query optimization for project list with sample and task counts
        This is a critical performance test
        """
        # Create test projects with samples and tasks
        for i in range(5):
            project = Project(
                name=f"Project {i}",
                description=f"Description {i}",
                user_id=test_user.id
            )
            db_session.add(project)
            db_session.commit()
            db_session.refresh(project)

            # Add samples to project
            for j in range(3):
                sample = Sample(
                    sample_id=f"Sample {i}-{j}",
                    project_id=project.id
                )
                db_session.add(sample)

            # Add tasks to project
            for k in range(2):
                task = PipelineTask(
                    task_name=f"Task {i}-{k}",
                    project_id=project.id,
                    status="pending"
                )
                db_session.add(task)

        db_session.commit()

        # Get project list
        response = client.get("/api/v1/projects", headers=auth_headers)
        assert response.status_code == 200
        data = response.json()

        # Verify response format
        assert data["total"] == 5
        assert len(data["items"]) == 5

        # Verify counts are included (N+1 optimization test)
        for project in data["items"]:
            assert "sample_count" in project
            assert "task_count" in project
            assert project["sample_count"] == 3
            assert project["task_count"] == 2

    def test_get_project_detail(self, client, auth_headers, db_session, test_user):
        """Test getting project details"""
        # Create a project
        project = Project(
            name="Detail Test Project",
            description="Test Description",
            user_id=test_user.id
        )
        db_session.add(project)
        db_session.commit()
        db_session.refresh(project)

        # Get project details
        response = client.get(f"/api/v1/projects/{project.id}", headers=auth_headers)
        assert response.status_code == 200
        data = response.json()
        assert data["id"] == str(project.id)
        assert data["name"] == "Detail Test Project"

    def test_update_project(self, client, auth_headers, db_session, test_user):
        """Test updating a project"""
        # Create a project
        project = Project(
            name="Original Name",
            user_id=test_user.id
        )
        db_session.add(project)
        db_session.commit()
        db_session.refresh(project)

        # Update project
        response = client.put(
            f"/api/v1/projects/{project.id}",
            headers=auth_headers,
            json={"name": "Updated Name", "description": "New Description"}
        )
        assert response.status_code == 200
        data = response.json()
        assert data["name"] == "Updated Name"
        assert data["description"] == "New Description"

    def test_delete_project(self, client, auth_headers, db_session, test_user):
        """Test deleting a project"""
        # Create a project
        project = Project(
            name="To Delete",
            user_id=test_user.id
        )
        db_session.add(project)
        db_session.commit()
        project_id = project.id

        # Delete project
        response = client.delete(f"/api/v1/projects/{project_id}", headers=auth_headers)
        assert response.status_code in [200, 204]

        # Verify deletion
        response = client.get(f"/api/v1/projects/{project_id}", headers=auth_headers)
        assert response.status_code == 404

    def test_project_permission_denied(self, client, auth_headers, db_session):
        """Test that users cannot access other users' projects"""
        # Create another user and their project
        from app.models.user import User
        from app.core.security import get_password_hash

        other_user = User(
            username="otheruser",
            email="other@example.com",
            password_hash=get_password_hash("password"),
            is_active=True
        )
        db_session.add(other_user)
        db_session.commit()
        db_session.refresh(other_user)

        project = Project(
            name="Other User's Project",
            user_id=other_user.id
        )
        db_session.add(project)
        db_session.commit()

        # Try to access with test_user credentials
        response = client.get(f"/api/v1/projects/{project.id}", headers=auth_headers)
        assert response.status_code == 404  # Not found (access denied)

    def test_get_project_stats(self, client, auth_headers, db_session, test_user):
        """Test getting project statistics"""
        # Create some test projects
        for status in ["active", "completed", "archived"]:
            project = Project(
                name=f"Project {status}",
                status=status,
                user_id=test_user.id
            )
            db_session.add(project)
        db_session.commit()

        response = client.get("/api/v1/projects/stats", headers=auth_headers)
        assert response.status_code == 200
        data = response.json()
        assert "total_projects" in data
        assert data["total_projects"] >= 3


class TestProjectsAPIErrorHandling:
    """Test error handling in Projects API"""

    def test_create_project_missing_required_fields(self, client, auth_headers):
        """Test validation error when required fields are missing"""
        response = client.post(
            "/api/v1/projects",
            headers=auth_headers,
            json={}  # Missing required 'name' field
        )
        assert response.status_code == 422  # Validation Error
        data = response.json()
        assert "error" in data or "detail" in data

    def test_get_nonexistent_project(self, client, auth_headers):
        """Test getting a project that doesn't exist"""
        fake_uuid = "00000000-0000-0000-0000-000000000000"
        response = client.get(f"/api/v1/projects/{fake_uuid}", headers=auth_headers)
        assert response.status_code == 404
        data = response.json()
        assert "error" in data or "detail" in data

    def test_unauthorized_access(self, client):
        """Test accessing projects without authentication"""
        response = client.get("/api/v1/projects")
        assert response.status_code == 401  # Unauthorized
