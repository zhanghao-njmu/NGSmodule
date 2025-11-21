"""
File storage service for MinIO/S3
"""
import hashlib
import os
from pathlib import Path
from typing import BinaryIO, Optional
import aiofiles
from minio import Minio
from minio.error import S3Error

from app.core.config import settings


class StorageService:
    """File storage service using MinIO/S3"""

    def __init__(self):
        self.client = Minio(
            settings.MINIO_URL,
            access_key=settings.MINIO_ACCESS_KEY,
            secret_key=settings.MINIO_SECRET_KEY,
            secure=settings.MINIO_SECURE
        )
        self.bucket_name = settings.MINIO_BUCKET
        self._ensure_bucket()

    def _ensure_bucket(self):
        """Ensure bucket exists"""
        try:
            if not self.client.bucket_exists(self.bucket_name):
                self.client.make_bucket(self.bucket_name)
        except S3Error as e:
            print(f"Error ensuring bucket: {e}")

    async def save_upload_file(
        self,
        file: BinaryIO,
        user_id: str,
        filename: str,
        project_id: Optional[str] = None
    ) -> str:
        """
        Save uploaded file to MinIO

        Args:
            file: File object
            user_id: User ID
            filename: Original filename
            project_id: Optional project ID

        Returns:
            Object path in storage
        """
        # Create object path
        if project_id:
            object_path = f"{user_id}/{project_id}/{filename}"
        else:
            object_path = f"{user_id}/uploads/{filename}"

        # Upload to MinIO
        try:
            # Get file size
            file.seek(0, os.SEEK_END)
            file_size = file.tell()
            file.seek(0)

            # Upload
            self.client.put_object(
                self.bucket_name,
                object_path,
                file,
                file_size
            )

            return object_path

        except S3Error as e:
            raise Exception(f"Error uploading file: {str(e)}")

    async def get_file(self, object_path: str) -> BinaryIO:
        """
        Get file from storage

        Args:
            object_path: Path to object in storage

        Returns:
            File stream
        """
        try:
            response = self.client.get_object(self.bucket_name, object_path)
            return response
        except S3Error as e:
            raise Exception(f"Error retrieving file: {str(e)}")

    async def delete_file(self, object_path: str):
        """
        Delete file from storage

        Args:
            object_path: Path to object in storage
        """
        try:
            self.client.remove_object(self.bucket_name, object_path)
        except S3Error as e:
            raise Exception(f"Error deleting file: {str(e)}")

    async def calculate_md5(self, file_path: str) -> str:
        """
        Calculate MD5 checksum of a file

        Args:
            file_path: Path to file

        Returns:
            MD5 hex digest
        """
        md5_hash = hashlib.md5()
        async with aiofiles.open(file_path, "rb") as f:
            while chunk := await f.read(8192):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()

    def get_presigned_url(
        self,
        object_path: str,
        expires_minutes: int = 60
    ) -> str:
        """
        Get presigned URL for direct download

        Args:
            object_path: Path to object
            expires_minutes: URL expiration time in minutes

        Returns:
            Presigned URL
        """
        try:
            from datetime import timedelta
            url = self.client.presigned_get_object(
                self.bucket_name,
                object_path,
                expires=timedelta(minutes=expires_minutes)
            )
            return url
        except S3Error as e:
            raise Exception(f"Error generating presigned URL: {str(e)}")


# Singleton instance
storage_service = StorageService()
