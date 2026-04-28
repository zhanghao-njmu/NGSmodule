"""
File Service - Business logic for file management
"""

import hashlib
import io
import os
from typing import List, Optional
from uuid import UUID

from fastapi import HTTPException, UploadFile, status
from sqlalchemy.orm import Session

from app.core.config import settings
from app.models.file import File
from app.models.project import Project
from app.models.sample import Sample
from app.models.user import User
from app.services.storage import storage_service


class FileService:
    """Service class for file-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize FileService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(self, file_id: UUID, user_id: UUID) -> Optional[File]:
        """
        Get file by ID for a specific user

        Args:
            file_id: File UUID
            user_id: User UUID (for authorization)

        Returns:
            File if found and belongs to user, None otherwise
        """
        file = (
            self.db.query(File)
            .join(Sample)
            .join(Project)
            .filter(File.id == file_id, Project.user_id == user_id)
            .first()
        )

        return file

    def get_by_id_or_raise(self, file_id: UUID, user_id: UUID) -> File:
        """
        Get file by ID or raise 404 error

        Args:
            file_id: File UUID
            user_id: User UUID (for authorization)

        Returns:
            File

        Raises:
            HTTPException: If file not found or unauthorized
        """
        file = self.get_by_id(file_id, user_id)
        if not file:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"File {file_id} not found")
        return file

    def list_files(
        self,
        user_id: UUID,
        sample_id: Optional[UUID] = None,
        project_id: Optional[UUID] = None,
        file_type: Optional[str] = None,
    ) -> List[File]:
        """
        List files for a user with optional filters

        Args:
            user_id: User UUID
            sample_id: Optional sample filter
            project_id: Optional project filter
            file_type: Optional file type filter

        Returns:
            List of files
        """
        query = self.db.query(File).join(Sample).join(Project).filter(Project.user_id == user_id)

        # Apply filters
        if sample_id:
            query = query.filter(File.sample_id == sample_id)

        if project_id:
            query = query.filter(Project.id == project_id)

        if file_type:
            query = query.filter(File.file_type == file_type)

        files = query.order_by(File.created_at.desc()).all()

        return files

    # ============= CREATE OPERATIONS =============

    async def upload(self, user_id: UUID, sample_id: UUID, file: UploadFile) -> File:
        """
        Upload a file and associate with sample

        Args:
            user_id: User UUID (for authorization)
            sample_id: Sample UUID to associate file with
            file: File to upload

        Returns:
            Created file record

        Raises:
            HTTPException: If sample not found, file type not allowed, or storage quota exceeded
        """
        # Verify sample belongs to user
        sample = self.db.query(Sample).join(Project).filter(Sample.id == sample_id, Project.user_id == user_id).first()

        if not sample:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Sample not found")

        # Get user
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="User not found")

        # Check file extension
        file_ext = os.path.splitext(file.filename)[1].lower()
        allowed_extensions = settings.ALLOWED_EXTENSIONS

        if file_ext not in allowed_extensions:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"File type '{file_ext}' not allowed. Allowed: {', '.join(allowed_extensions)}",
            )

        # Check storage quota
        file_size = 0
        try:
            file_contents = await file.read()
            file_size = len(file_contents)

            if user.storage_used + file_size > user.storage_quota:
                raise HTTPException(
                    status_code=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE, detail="Storage quota exceeded"
                )

            # Reset file pointer
            await file.seek(0)

        except HTTPException:
            raise
        except Exception as e:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Error reading file: {str(e)}")

        # Save to storage
        try:
            file_stream = io.BytesIO(file_contents)
            object_path = await storage_service.save_upload_file(
                file=file_stream, user_id=str(user_id), filename=file.filename, project_id=str(sample.project_id)
            )

            # Calculate MD5
            md5_hash = hashlib.md5(file_contents).hexdigest()

            # Create file record
            file_record = File(
                sample_id=sample_id,
                filename=file.filename,
                file_path=object_path,
                file_type=file_ext.lstrip("."),
                file_size=file_size,
                md5_checksum=md5_hash,
                upload_status="completed",
            )

            self.db.add(file_record)

            # Update user storage
            user.storage_used += file_size

            self.db.commit()
            self.db.refresh(file_record)

            return file_record

        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Error uploading file: {str(e)}"
            )

    # ============= DOWNLOAD OPERATIONS =============

    def get_download_url(self, file_id: UUID, user_id: UUID) -> str:
        """
        Get presigned download URL for a file

        Args:
            file_id: File UUID
            user_id: User UUID (for authorization)

        Returns:
            Presigned download URL

        Raises:
            HTTPException: If file not found or error generating URL
        """
        file = self.get_by_id_or_raise(file_id, user_id)

        try:
            # Get presigned URL (redirect to MinIO)
            presigned_url = storage_service.get_presigned_url(file.file_path)
            return presigned_url

        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Error generating download URL: {str(e)}"
            )

    # ============= DELETE OPERATIONS =============

    async def delete(self, file_id: UUID, user_id: UUID) -> None:
        """
        Delete file from storage and database

        Args:
            file_id: File UUID
            user_id: User UUID (for authorization)

        Raises:
            HTTPException: If file not found or error deleting
        """
        file = self.get_by_id_or_raise(file_id, user_id)

        # Get user
        user = self.db.query(User).filter(User.id == user_id).first()
        if not user:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="User not found")

        try:
            # Delete from storage
            await storage_service.delete_file(file.file_path)

            # Update user storage
            if file.file_size:
                user.storage_used = max(0, user.storage_used - file.file_size)

            # Delete record
            self.db.delete(file)
            self.db.commit()

        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Error deleting file: {str(e)}"
            )

    # ============= HELPER METHODS =============

    def verify_file_checksum(self, file_id: UUID, user_id: UUID, provided_checksum: str) -> bool:
        """
        Verify file MD5 checksum

        Args:
            file_id: File UUID
            user_id: User UUID (for authorization)
            provided_checksum: Checksum to verify

        Returns:
            True if checksums match, False otherwise

        Raises:
            HTTPException: If file not found
        """
        file = self.get_by_id_or_raise(file_id, user_id)

        if not file.md5_checksum:
            return False

        return file.md5_checksum == provided_checksum

    def get_files_by_sample(self, sample_id: UUID, user_id: UUID) -> List[File]:
        """
        Get all files for a specific sample

        Args:
            sample_id: Sample UUID
            user_id: User UUID (for authorization)

        Returns:
            List of files

        Raises:
            HTTPException: If sample not found
        """
        # Verify sample belongs to user
        sample = self.db.query(Sample).join(Project).filter(Sample.id == sample_id, Project.user_id == user_id).first()

        if not sample:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Sample not found")

        files = self.db.query(File).filter(File.sample_id == sample_id).order_by(File.created_at.desc()).all()

        return files
