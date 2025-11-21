"""
Celery tasks for NGS pipeline execution
"""
import subprocess
import os
from pathlib import Path
from typing import Dict, Any
from uuid import UUID
import asyncio

from app.workers.celery_app import celery_app
from app.core.config import settings
from app.core.database import SessionLocal
from app.models.task import PipelineTask


def send_websocket_update(task_id: str, status: str, progress: float, message: str = ""):
    """
    Send WebSocket update for task progress

    Note: This runs in a new event loop since Celery workers don't have an async loop
    """
    try:
        from app.api.v1.websocket import manager

        # Create new event loop for this thread
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        try:
            loop.run_until_complete(
                manager.send_task_status(task_id, status, progress, message)
            )
        finally:
            loop.close()
    except Exception as e:
        print(f"Error sending WebSocket update: {e}")


@celery_app.task(bind=True)
def run_ngs_pipeline(
    self,
    task_id: str,
    pipeline_script: str,
    config: Dict[str, Any]
):
    """
    Run NGS pipeline script

    Args:
        self: Celery task instance
        task_id: PipelineTask ID
        pipeline_script: Path to pipeline script
        config: Pipeline configuration

    Returns:
        Dict with task results
    """
    db = SessionLocal()

    try:
        # Get task
        task = db.query(PipelineTask).filter(PipelineTask.id == task_id).first()
        if not task:
            raise Exception(f"Task {task_id} not found")

        # Update task status
        task.status = "running"
        task.celery_task_id = self.request.id
        db.commit()

        # Update progress
        self.update_state(
            state="PROGRESS",
            meta={"progress": 0, "status": "Starting pipeline..."}
        )

        # Send WebSocket update
        send_websocket_update(task_id, "running", 0.0, "Starting pipeline...")

        # Prepare command
        script_path = Path(settings.NGS_PIPELINE_DIR) / pipeline_script
        if not script_path.exists():
            raise Exception(f"Pipeline script not found: {script_path}")

        # Create log file
        log_dir = Path(settings.NGS_WORK_DIR) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / f"{task_id}.log"
        task.log_file_path = str(log_file)
        db.commit()

        # Build command
        cmd = ["bash", str(script_path)]

        # Add configuration parameters
        for key, value in config.items():
            cmd.extend([f"--{key}", str(value)])

        # Run pipeline
        with open(log_file, "w") as log:
            process = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                cwd=settings.NGS_PIPELINE_DIR,
                env={**os.environ, "NGS_TASK_ID": task_id}
            )

            # Wait for completion
            return_code = process.wait()

            if return_code != 0:
                raise Exception(f"Pipeline failed with return code {return_code}")

        # Update task status
        task.status = "completed"
        task.progress = 100.0
        db.commit()

        self.update_state(
            state="SUCCESS",
            meta={"progress": 100, "status": "Pipeline completed successfully"}
        )

        # Send WebSocket update
        send_websocket_update(task_id, "completed", 100.0, "Pipeline completed successfully")

        return {
            "task_id": task_id,
            "status": "completed",
            "log_file": str(log_file)
        }

    except Exception as e:
        # Update task status
        if task:
            task.status = "failed"
            task.error_message = str(e)
            db.commit()

            # Send WebSocket update
            send_websocket_update(task_id, "failed", task.progress, f"Pipeline failed: {str(e)}")

        self.update_state(
            state="FAILURE",
            meta={"error": str(e)}
        )

        raise

    finally:
        db.close()


@celery_app.task
def cleanup_old_files(days: int = 30):
    """
    Clean up old temporary files

    Args:
        days: Delete files older than this many days
    """
    import time
    from datetime import datetime, timedelta

    work_dir = Path(settings.NGS_WORK_DIR)
    cutoff_time = time.time() - (days * 86400)

    deleted_count = 0

    for file_path in work_dir.rglob("*"):
        if file_path.is_file():
            if file_path.stat().st_mtime < cutoff_time:
                try:
                    file_path.unlink()
                    deleted_count += 1
                except Exception as e:
                    print(f"Error deleting {file_path}: {e}")

    return {
        "deleted_files": deleted_count,
        "cutoff_date": datetime.fromtimestamp(cutoff_time).isoformat()
    }
