"""
Project Service - Business logic for project management
"""
from typing import List, Optional, Dict, Any, Tuple
from uuid import UUID
from sqlalchemy.orm import Session
from sqlalchemy import func, or_
from fastapi import HTTPException, status

from app.models.project import Project
from app.models.sample import Sample
from app.models.task import PipelineTask
from app.models.user import User
from app.schemas.project import (
    ProjectCreate,
    ProjectUpdate,
    ProjectStats,
)


class ProjectService:
    """Service class for project-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize ProjectService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(
        self,
        project_id: UUID,
        user_id: UUID,
        is_admin: bool = False,
    ) -> Optional[Project]:
        """
        Get project by ID.

        Args:
            project_id: Project UUID
            user_id: Caller's user UUID (used for ownership check unless admin)
            is_admin: If True, bypass the ownership check

        Returns:
            Project if found (and accessible), None otherwise

        Note:
            Soft-deleted projects (status='deleted') are excluded for both
            owners and admins.
        """
        query = self.db.query(Project).filter(
            Project.id == project_id,
            Project.status != "deleted",
        )
        if not is_admin:
            query = query.filter(Project.user_id == user_id)

        project = query.first()

        if project:
            # Add computed fields
            project.sample_count = self._get_sample_count(project_id)
            project.task_count = self._get_task_count(project_id)

        return project

    def get_by_id_or_raise(
        self,
        project_id: UUID,
        user_id: UUID,
        is_admin: bool = False,
    ) -> Project:
        """
        Get project by ID or raise 404 error

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)
            is_admin: If True, bypass ownership check

        Returns:
            Project

        Raises:
            HTTPException: If project not found or unauthorized
        """
        project = self.get_by_id(project_id, user_id, is_admin=is_admin)
        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Project {project_id} not found"
            )
        return project

    def list_projects(
        self,
        user_id: UUID,
        skip: int = 0,
        limit: int = 20,
        status_filter: Optional[str] = None,
        project_type: Optional[str] = None,
        search: Optional[str] = None,
    ) -> Tuple[List[Project], int]:
        """
        List projects for a user with pagination and filters

        Args:
            user_id: User UUID
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            status_filter: Optional status filter
            project_type: Optional project type filter
            search: Optional search query (name)

        Returns:
            Tuple of (projects list, total count)
        """
        query = self.db.query(Project).filter(Project.user_id == user_id)

        # Apply filters
        if status_filter:
            query = query.filter(Project.status == status_filter)
        if project_type:
            query = query.filter(Project.project_type == project_type)
        if search:
            query = query.filter(Project.name.ilike(f"%{search}%"))

        # Get total count before pagination
        total = query.count()

        # Apply pagination and ordering
        projects = query.order_by(Project.created_at.desc()).offset(skip).limit(limit).all()

        # Add computed fields efficiently (avoid N+1 queries)
        if projects:
            self._add_computed_fields(projects)

        return projects, total

    def get_stats(self, user_id: UUID) -> ProjectStats:
        """
        Get project statistics for a user

        Args:
            user_id: User UUID

        Returns:
            ProjectStats with counts
        """
        # Count projects by status
        total_projects = self.db.query(Project).filter(
            Project.user_id == user_id
        ).count()

        active_projects = self.db.query(Project).filter(
            Project.user_id == user_id,
            Project.status == "active"
        ).count()

        archived_projects = self.db.query(Project).filter(
            Project.user_id == user_id,
            Project.status == "archived"
        ).count()

        # Count samples across all user projects
        total_samples = self.db.query(Sample).join(Project).filter(
            Project.user_id == user_id
        ).count()

        # Count tasks across all user projects
        total_tasks = self.db.query(PipelineTask).join(Project).filter(
            Project.user_id == user_id
        ).count()

        running_tasks = self.db.query(PipelineTask).join(Project).filter(
            Project.user_id == user_id,
            PipelineTask.status == "running"
        ).count()

        completed_tasks = self.db.query(PipelineTask).join(Project).filter(
            Project.user_id == user_id,
            PipelineTask.status == "completed"
        ).count()

        failed_tasks = self.db.query(PipelineTask).join(Project).filter(
            Project.user_id == user_id,
            PipelineTask.status == "failed"
        ).count()

        return ProjectStats(
            total_projects=total_projects,
            active_projects=active_projects,
            archived_projects=archived_projects,
            total_samples=total_samples,
            total_tasks=total_tasks,
            running_tasks=running_tasks,
            completed_tasks=completed_tasks,
            failed_tasks=failed_tasks,
        )

    # ============= CREATE OPERATIONS =============

    def create(self, user_id: UUID, project_data: ProjectCreate) -> Project:
        """
        Create a new project

        Args:
            user_id: User UUID
            project_data: Project creation data

        Returns:
            Created project

        Raises:
            HTTPException: If project name already exists for user
        """
        # Check for duplicate name
        existing = self.db.query(Project).filter(
            Project.user_id == user_id,
            Project.name == project_data.name
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Project '{project_data.name}' already exists"
            )

        # Create project
        project = Project(
            user_id=user_id,
            name=project_data.name,
            description=project_data.description,
            project_type=project_data.project_type,
            config=project_data.config,
            status="active"
        )

        self.db.add(project)
        self.db.commit()
        self.db.refresh(project)

        # Add computed fields
        project.sample_count = 0
        project.task_count = 0

        return project

    # ============= UPDATE OPERATIONS =============

    def update(
        self,
        project_id: UUID,
        user_id: UUID,
        update_data: ProjectUpdate
    ) -> Project:
        """
        Update an existing project

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)
            update_data: Project update data

        Returns:
            Updated project

        Raises:
            HTTPException: If project not found or name conflict
        """
        # Get project (raises 404 if not found)
        project = self.get_by_id_or_raise(project_id, user_id)

        # Check for name conflict if name is being changed
        if update_data.name and update_data.name != project.name:
            existing = self.db.query(Project).filter(
                Project.user_id == user_id,
                Project.name == update_data.name,
                Project.id != project_id
            ).first()

            if existing:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"Project '{update_data.name}' already exists"
                )

        # Update fields
        update_dict = update_data.model_dump(exclude_unset=True)
        for field, value in update_dict.items():
            setattr(project, field, value)

        self.db.commit()
        self.db.refresh(project)

        # Add computed fields
        project.sample_count = self._get_sample_count(project_id)
        project.task_count = self._get_task_count(project_id)

        return project

    def archive(self, project_id: UUID, user_id: UUID) -> Project:
        """
        Archive a project

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)

        Returns:
            Archived project
        """
        project = self.get_by_id_or_raise(project_id, user_id)
        project.status = "archived"
        self.db.commit()
        self.db.refresh(project)
        return project

    def restore(self, project_id: UUID, user_id: UUID) -> Project:
        """
        Restore an archived project

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)

        Returns:
            Restored project
        """
        project = self.get_by_id_or_raise(project_id, user_id)
        project.status = "active"
        self.db.commit()
        self.db.refresh(project)
        return project

    # ============= DELETE OPERATIONS =============

    def delete(self, project_id: UUID, user_id: UUID) -> None:
        """
        Delete a project

        Args:
            project_id: Project UUID
            user_id: User UUID (for authorization)

        Raises:
            HTTPException: If project not found or has dependencies
        """
        project = self.get_by_id_or_raise(project_id, user_id)

        # Check for dependencies (samples, tasks)
        sample_count = self._get_sample_count(project_id)
        task_count = self._get_task_count(project_id)

        if sample_count > 0 or task_count > 0:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Cannot delete project with {sample_count} samples and {task_count} tasks. "
                       "Delete or move them first, or use archive instead."
            )

        self.db.delete(project)
        self.db.commit()

    # ============= HELPER METHODS =============

    def _get_sample_count(self, project_id: UUID) -> int:
        """Get sample count for a project"""
        return self.db.query(Sample).filter(Sample.project_id == project_id).count()

    def _get_task_count(self, project_id: UUID) -> int:
        """Get task count for a project"""
        return self.db.query(PipelineTask).filter(PipelineTask.project_id == project_id).count()

    def _add_computed_fields(self, projects: List[Project]) -> None:
        """
        Add computed fields to a list of projects efficiently (avoids N+1 queries)

        Args:
            projects: List of Project objects
        """
        project_ids = [p.id for p in projects]

        # Get sample counts in one query
        sample_counts = dict(
            self.db.query(Sample.project_id, func.count(Sample.id))
            .filter(Sample.project_id.in_(project_ids))
            .group_by(Sample.project_id)
            .all()
        )

        # Get task counts in one query
        task_counts = dict(
            self.db.query(PipelineTask.project_id, func.count(PipelineTask.id))
            .filter(PipelineTask.project_id.in_(project_ids))
            .group_by(PipelineTask.project_id)
            .all()
        )

        # Assign counts to projects
        for project in projects:
            project.sample_count = sample_counts.get(project.id, 0)
            project.task_count = task_counts.get(project.id, 0)
