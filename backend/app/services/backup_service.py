"""
Backup Service
Real implementation using pg_dump for database backups and tar for files.
"""

import hashlib
import logging
import os
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

from sqlalchemy.orm import Session

from app.core.config import settings
from app.core.datetime_utils import utc_now_naive
from app.models.backup import SystemBackup

logger = logging.getLogger(__name__)


class BackupService:
    """Service for system backup operations"""

    def __init__(self, db: Session):
        self.db = db
        self.backup_dir = Path(settings.BACKUP_DIR)
        self.backup_dir.mkdir(parents=True, exist_ok=True)

    def _parse_database_url(self) -> dict:
        """Parse DATABASE_URL into pg_dump parameters"""
        parsed = urlparse(settings.DATABASE_URL)
        return {
            "host": parsed.hostname or "localhost",
            "port": str(parsed.port or 5432),
            "user": parsed.username or "postgres",
            "password": parsed.password or "",
            "database": parsed.path.lstrip("/") if parsed.path else "ngsmodule",
        }

    def _calculate_checksum(self, file_path: Path) -> str:
        """Calculate SHA256 checksum of a file"""
        sha256 = hashlib.sha256()
        try:
            with open(file_path, "rb") as f:
                for chunk in iter(lambda: f.read(65536), b""):
                    sha256.update(chunk)
            return sha256.hexdigest()
        except Exception as e:
            logger.error(f"Failed to calculate checksum for {file_path}: {e}")
            return ""

    def _backup_database(self, output_path: Path, compress: bool) -> bool:
        """
        Backup database using pg_dump

        Returns True on success, False on failure.
        """
        db_params = self._parse_database_url()

        env = os.environ.copy()
        env["PGPASSWORD"] = db_params["password"]

        cmd = [
            "pg_dump",
            "-h",
            db_params["host"],
            "-p",
            db_params["port"],
            "-U",
            db_params["user"],
            "-d",
            db_params["database"],
            "-F",
            "c" if compress else "p",  # custom (compressed) or plain format
            "-f",
            str(output_path),
        ]

        try:
            result = subprocess.run(
                cmd,
                env=env,
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour max
            )
            if result.returncode != 0:
                logger.error(f"pg_dump failed: {result.stderr}")
                return False
            return True
        except subprocess.TimeoutExpired:
            logger.error("pg_dump timed out")
            return False
        except FileNotFoundError:
            logger.error("pg_dump not found. Is PostgreSQL client installed?")
            return False
        except Exception as e:
            logger.error(f"pg_dump error: {e}")
            return False

    def _backup_files(self, output_path: Path, compress: bool) -> bool:
        """
        Backup user files using tar

        Returns True on success, False on failure.
        """
        upload_dir = Path(settings.UPLOAD_DIR)
        if not upload_dir.exists():
            logger.warning(f"Upload directory {upload_dir} does not exist, skipping file backup")
            # Create empty archive
            output_path.touch()
            return True

        cmd = ["tar"]
        if compress:
            cmd.extend(["-czf", str(output_path)])
        else:
            cmd.extend(["-cf", str(output_path)])
        cmd.extend(["-C", str(upload_dir.parent), upload_dir.name])

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200,  # 2 hours max
            )
            if result.returncode != 0:
                logger.error(f"tar failed: {result.stderr}")
                return False
            return True
        except subprocess.TimeoutExpired:
            logger.error("tar timed out")
            return False
        except Exception as e:
            logger.error(f"tar error: {e}")
            return False

    def create_backup(
        self,
        backup_type: str,
        admin_user_id: str,
        description: Optional[str] = None,
        compress: bool = True,
    ) -> SystemBackup:
        """
        Create a system backup synchronously

        For production use, this should be called from a Celery task.
        """
        timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        extension = "dump" if compress else "sql"

        # Determine file name based on backup type
        if backup_type == "database_only":
            file_name = f"db_backup_{timestamp}.{extension}"
        elif backup_type == "files_only":
            file_name = f"files_backup_{timestamp}.tar{'.gz' if compress else ''}"
        else:
            file_name = f"full_backup_{timestamp}.tar{'.gz' if compress else ''}"

        file_path = self.backup_dir / file_name

        # Create backup record
        backup = SystemBackup(
            backup_type=backup_type,
            status="in_progress",
            file_path=str(file_path),
            file_name=file_name,
            compressed=compress,
            description=description,
            created_by=admin_user_id,
            started_at=utc_now_naive(),
            expires_at=utc_now_naive() + timedelta(days=settings.BACKUP_RETENTION_DAYS),
        )
        self.db.add(backup)
        self.db.commit()
        self.db.refresh(backup)

        try:
            success = False

            if backup_type == "database_only":
                success = self._backup_database(file_path, compress)

            elif backup_type == "files_only":
                success = self._backup_files(file_path, compress)

            elif backup_type in ("full", "incremental"):
                # Full backup: database + files
                db_path = self.backup_dir / f"db_{timestamp}.dump"
                files_path = self.backup_dir / f"files_{timestamp}.tar.gz"

                db_success = self._backup_database(db_path, True)
                files_success = self._backup_files(files_path, True)

                if db_success and files_success:
                    # Combine into single archive
                    cmd = [
                        "tar",
                        "-czf" if compress else "-cf",
                        str(file_path),
                        "-C",
                        str(self.backup_dir),
                        db_path.name,
                        files_path.name,
                    ]
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    success = result.returncode == 0

                    # Cleanup intermediate files
                    db_path.unlink(missing_ok=True)
                    files_path.unlink(missing_ok=True)
                else:
                    success = False

            if success and file_path.exists():
                backup.status = "completed"
                backup.size = file_path.stat().st_size
                backup.checksum = self._calculate_checksum(file_path)
                backup.completed_at = utc_now_naive()
            else:
                backup.status = "failed"
                backup.error_message = "Backup operation failed. Check logs for details."

        except Exception as e:
            logger.error(f"Backup failed: {e}", exc_info=True)
            backup.status = "failed"
            backup.error_message = str(e)

        self.db.commit()
        self.db.refresh(backup)
        return backup

    def list_backups(self, limit: int = 100) -> list:
        """List all backups, newest first"""
        return self.db.query(SystemBackup).order_by(SystemBackup.created_at.desc()).limit(limit).all()

    def get_backup(self, backup_id: str) -> Optional[SystemBackup]:
        """Get backup by ID"""
        return self.db.query(SystemBackup).filter(SystemBackup.id == backup_id).first()

    def delete_backup(self, backup_id: str) -> bool:
        """Delete a backup file and database record"""
        backup = self.get_backup(backup_id)
        if not backup:
            return False

        # Delete file
        if backup.file_path:
            try:
                Path(backup.file_path).unlink(missing_ok=True)
            except Exception as e:
                logger.error(f"Failed to delete backup file {backup.file_path}: {e}")

        # Delete record
        self.db.delete(backup)
        self.db.commit()
        return True

    def cleanup_expired_backups(self) -> int:
        """Delete backups past their expiration date"""
        now = utc_now_naive()
        expired = self.db.query(SystemBackup).filter(SystemBackup.expires_at < now).all()

        count = 0
        for backup in expired:
            if self.delete_backup(str(backup.id)):
                count += 1

        return count
