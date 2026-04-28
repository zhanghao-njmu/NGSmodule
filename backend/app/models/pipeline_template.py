"""
Pipeline Template Model

Represents predefined pipeline templates for common NGS analysis workflows.
"""

import uuid

from sqlalchemy import JSON, Boolean, Column, String, Text

from app.core.database import Base
from app.core.types import UUID


class PipelineTemplate(Base):
    """Pipeline template model for predefined workflows"""

    __tablename__ = "pipeline_templates"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    name = Column(String(100), nullable=False, unique=True, index=True)
    display_name = Column(String(200), nullable=False)
    description = Column(Text)
    category = Column(String(50), nullable=False, index=True)  # 'RNA-seq', 'DNA-seq', 'scRNA-seq', etc.

    # Pipeline script information
    script_name = Column(String(200), nullable=False)  # e.g., 'preAlignmentQC', 'Alignment'
    script_path = Column(String(500))  # Relative path to the script

    # Default parameters (JSON)
    default_params = Column(JSON, default={})

    # Parameter schema (defines what parameters are required/optional)
    param_schema = Column(JSON, default={})

    # Execution requirements
    estimated_time = Column(String(50))  # e.g., '2-4 hours'
    min_memory_gb = Column(String(20))  # e.g., '16GB'
    min_cpu_cores = Column(String(20))  # e.g., '4'

    # Status
    is_active = Column(Boolean, default=True)
    is_builtin = Column(Boolean, default=True)  # Built-in templates cannot be deleted

    # Ordering and grouping
    sort_order = Column(String(10), default="0")
    tags = Column(JSON, default=[])  # e.g., ['quality-control', 'alignment']

    def __repr__(self):
        return f"<PipelineTemplate {self.name} - {self.display_name}>"
