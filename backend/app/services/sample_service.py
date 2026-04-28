"""
Sample Service - Business logic for sample management
"""

import csv
import io
from typing import List, Optional, Tuple
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy import func
from sqlalchemy.orm import Session

from app.models.file import File as FileModel
from app.models.project import Project
from app.models.sample import Sample
from app.schemas.sample import (
    SampleBase,
    SampleCreate,
    SampleUpdate,
)


class SampleService:
    """Service class for sample-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize SampleService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(self, sample_id: UUID, user_id: UUID) -> Optional[Sample]:
        """
        Get sample by ID for a specific user

        Args:
            sample_id: Sample UUID
            user_id: User UUID (for authorization)

        Returns:
            Sample if found and belongs to user, None otherwise
        """
        sample = self.db.query(Sample).join(Project).filter(Sample.id == sample_id, Project.user_id == user_id).first()

        if sample:
            # Add computed fields
            sample.file_count = self._get_file_count(sample_id)

        return sample

    def get_by_id_or_raise(self, sample_id: UUID, user_id: UUID) -> Sample:
        """
        Get sample by ID or raise 404 error

        Args:
            sample_id: Sample UUID
            user_id: User UUID (for authorization)

        Returns:
            Sample

        Raises:
            HTTPException: If sample not found or unauthorized
        """
        sample = self.get_by_id(sample_id, user_id)
        if not sample:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Sample {sample_id} not found")
        return sample

    def list_samples(
        self,
        user_id: UUID,
        skip: int = 0,
        limit: int = 50,
        project_id: Optional[UUID] = None,
        group_name: Optional[str] = None,
    ) -> Tuple[List[Sample], int]:
        """
        List samples for a user with pagination and filters

        Args:
            user_id: User UUID
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            project_id: Optional project filter
            group_name: Optional group name filter

        Returns:
            Tuple of (list of samples, total count)
        """
        query = self.db.query(Sample).join(Project).filter(Project.user_id == user_id)

        # Apply filters
        if project_id:
            # Verify project belongs to user
            project = self.db.query(Project).filter(Project.id == project_id, Project.user_id == user_id).first()
            if not project:
                raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")
            query = query.filter(Sample.project_id == project_id)

        if group_name:
            query = query.filter(Sample.group_name == group_name)

        # Get total count
        total = query.count()

        # Get paginated results
        samples = query.order_by(Sample.created_at.desc()).offset(skip).limit(limit).all()

        # Add computed fields efficiently (prevent N+1 queries)
        if samples:
            self._add_computed_fields(samples)

        return samples, total

    # ============= CREATE OPERATIONS =============

    def create(self, user_id: UUID, sample_data: SampleCreate) -> Sample:
        """
        Create a new sample

        Args:
            user_id: User UUID (for authorization)
            sample_data: Sample creation data

        Returns:
            Created sample

        Raises:
            HTTPException: If project not found or sample_id already exists
        """
        # Verify project belongs to user
        project = (
            self.db.query(Project).filter(Project.id == sample_data.project_id, Project.user_id == user_id).first()
        )

        if not project:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

        # Check if sample_id already exists in this project
        existing = (
            self.db.query(Sample)
            .filter(Sample.project_id == sample_data.project_id, Sample.sample_id == sample_data.sample_id)
            .first()
        )

        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Sample '{sample_data.sample_id}' already exists in this project",
            )

        # Create sample
        sample = Sample(
            project_id=sample_data.project_id,
            sample_id=sample_data.sample_id,
            run_id=sample_data.run_id,
            group_name=sample_data.group_name,
            layout=sample_data.layout,
            batch_id=sample_data.batch_id,
            sample_metadata=sample_data.metadata or {},
        )

        self.db.add(sample)
        self.db.commit()
        self.db.refresh(sample)

        # Add computed fields
        sample.file_count = 0

        return sample

    def create_batch(self, user_id: UUID, project_id: UUID, samples_data: List[SampleBase]) -> int:
        """
        Create multiple samples at once

        Args:
            user_id: User UUID (for authorization)
            project_id: Project UUID
            samples_data: List of sample data

        Returns:
            Number of samples created

        Raises:
            HTTPException: If project not found or duplicate sample IDs
        """
        # Verify project belongs to user
        project = self.db.query(Project).filter(Project.id == project_id, Project.user_id == user_id).first()

        if not project:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

        # Get existing sample IDs
        existing_samples = self.db.query(Sample.sample_id).filter(Sample.project_id == project_id).all()
        existing_ids = {s.sample_id for s in existing_samples}

        # Check for duplicates
        new_ids = [s.sample_id for s in samples_data]
        duplicates = [sid for sid in new_ids if sid in existing_ids]

        if duplicates:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail=f"Duplicate sample IDs found: {', '.join(duplicates)}"
            )

        # Create samples
        samples_to_add = []
        for sample_data in samples_data:
            sample = Sample(
                project_id=project_id,
                sample_id=sample_data.sample_id,
                run_id=sample_data.run_id,
                group_name=sample_data.group_name,
                layout=sample_data.layout,
                batch_id=sample_data.batch_id,
                sample_metadata={},
            )
            samples_to_add.append(sample)

        self.db.add_all(samples_to_add)
        self.db.commit()

        return len(samples_to_add)

    def import_from_csv(self, user_id: UUID, project_id: UUID, csv_content: bytes) -> Tuple[int, str]:
        """
        Import samples from CSV file

        Args:
            user_id: User UUID (for authorization)
            project_id: Project UUID
            csv_content: CSV file content as bytes

        Returns:
            Tuple of (number of samples imported, project name)

        Raises:
            HTTPException: If project not found or CSV parsing errors
        """
        # Verify project belongs to user
        project = self.db.query(Project).filter(Project.id == project_id, Project.user_id == user_id).first()

        if not project:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Project not found")

        # Read CSV content
        try:
            decoded = csv_content.decode("utf-8")
            csv_reader = csv.DictReader(io.StringIO(decoded))

            # Validate headers
            required_headers = {"sample_id"}
            if not required_headers.issubset(csv_reader.fieldnames or []):
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"CSV must contain at least: {', '.join(required_headers)}",
                )

            # Get existing sample IDs
            existing_samples = self.db.query(Sample.sample_id).filter(Sample.project_id == project_id).all()
            existing_ids = {s.sample_id for s in existing_samples}

            # Parse and create samples
            samples_to_add = []
            for row in csv_reader:
                sample_id = row.get("sample_id", "").strip()
                if not sample_id:
                    continue

                if sample_id in existing_ids:
                    continue  # Skip duplicates

                sample = Sample(
                    project_id=project_id,
                    sample_id=sample_id,
                    run_id=row.get("run_id", "").strip() or None,
                    group_name=row.get("group_name", "").strip() or None,
                    layout=row.get("layout", "").strip() or None,
                    batch_id=row.get("batch_id", "").strip() or None,
                    sample_metadata={},
                )
                samples_to_add.append(sample)
                existing_ids.add(sample_id)

            if not samples_to_add:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail="No valid samples found in CSV or all samples already exist",
                )

            self.db.add_all(samples_to_add)
            self.db.commit()

            return len(samples_to_add), project.name

        except UnicodeDecodeError:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid CSV file encoding. Please use UTF-8 encoding."
            )
        except HTTPException:
            raise
        except Exception as e:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Error parsing CSV: {str(e)}")

    # ============= UPDATE OPERATIONS =============

    def update(self, sample_id: UUID, user_id: UUID, update_data: SampleUpdate) -> Sample:
        """
        Update sample information

        Args:
            sample_id: Sample UUID
            user_id: User UUID (for authorization)
            update_data: Update data

        Returns:
            Updated sample

        Raises:
            HTTPException: If sample not found or sample_id conflicts
        """
        sample = self.get_by_id_or_raise(sample_id, user_id)

        # Check if new sample_id conflicts
        if update_data.sample_id and update_data.sample_id != sample.sample_id:
            existing = (
                self.db.query(Sample)
                .filter(
                    Sample.project_id == sample.project_id,
                    Sample.sample_id == update_data.sample_id,
                    Sample.id != sample_id,
                )
                .first()
            )

            if existing:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail=f"Sample '{update_data.sample_id}' already exists in this project",
                )

        # Update fields
        update_dict = update_data.dict(exclude_unset=True)
        for field, value in update_dict.items():
            setattr(sample, field, value)

        self.db.commit()
        self.db.refresh(sample)

        # Add computed fields
        sample.file_count = self._get_file_count(sample.id)

        return sample

    # ============= DELETE OPERATIONS =============

    def delete(self, sample_id: UUID, user_id: UUID) -> None:
        """
        Delete sample and all associated files

        Args:
            sample_id: Sample UUID
            user_id: User UUID (for authorization)

        Raises:
            HTTPException: If sample not found
        """
        sample = self.get_by_id_or_raise(sample_id, user_id)

        # Delete sample (cascade will delete associated files)
        self.db.delete(sample)
        self.db.commit()

    # ============= HELPER METHODS =============

    def _get_file_count(self, sample_id: UUID) -> int:
        """
        Get file count for a sample

        Args:
            sample_id: Sample UUID

        Returns:
            File count
        """
        return self.db.query(FileModel).filter(FileModel.sample_id == sample_id).count()

    def _add_computed_fields(self, samples: List[Sample]) -> None:
        """
        Add computed fields efficiently (avoids N+1 queries)

        Args:
            samples: List of samples to add fields to
        """
        sample_ids = [s.id for s in samples]

        # Get file counts in ONE query (not N queries)
        file_counts = dict(
            self.db.query(FileModel.sample_id, func.count(FileModel.id))
            .filter(FileModel.sample_id.in_(sample_ids))
            .group_by(FileModel.sample_id)
            .all()
        )

        # Assign counts
        for sample in samples:
            sample.file_count = file_counts.get(sample.id, 0)
