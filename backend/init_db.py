"""
Database initialization script

This script initializes the database with all required tables.
Run this before starting the application for the first time.
"""
from app.core.database import Base, engine
from app.models.user import User
from app.models.project import Project
from app.models.sample import Sample
from app.models.file import File
from app.models.task import PipelineTask
from app.models.result import Result
from app.models.pipeline_template import PipelineTemplate


def init_db():
    """
    Create all database tables
    """
    print("Creating database tables...")
    Base.metadata.create_all(bind=engine)
    print("✅ Database tables created successfully!")


if __name__ == "__main__":
    init_db()
