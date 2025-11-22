"""
Pipeline Service - Business logic for pipeline template management
"""
from typing import List, Optional, Dict, Any, Tuple
from uuid import UUID
from sqlalchemy.orm import Session
from sqlalchemy import or_
from fastapi import HTTPException, status
from collections import defaultdict
from datetime import datetime

from app.models.pipeline_template import PipelineTemplate
from app.models.project import Project
from app.models.task import PipelineTask
from app.models.sample import Sample
from app.schemas.pipeline import (
    PipelineTemplateCreate,
    PipelineTemplateUpdate,
    PipelineTemplateCategory,
    PipelineExecuteRequest,
    PipelineBatchExecuteRequest,
    PipelineBatchExecuteResponse,
    ParameterRecommendationResponse,
)
from app.workers.pipeline_tasks import run_ngs_pipeline


class PipelineService:
    """Service class for pipeline template-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize PipelineService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(self, template_id: UUID) -> Optional[PipelineTemplate]:
        """
        Get pipeline template by ID

        Args:
            template_id: Template UUID

        Returns:
            Template if found, None otherwise
        """
        return self.db.query(PipelineTemplate).filter(
            PipelineTemplate.id == template_id
        ).first()

    def get_by_id_or_raise(self, template_id: UUID) -> PipelineTemplate:
        """
        Get pipeline template by ID or raise 404 error

        Args:
            template_id: Template UUID

        Returns:
            Template

        Raises:
            HTTPException: If template not found
        """
        template = self.get_by_id(template_id)
        if not template:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Pipeline template {template_id} not found"
            )
        return template

    def list_templates(
        self,
        category: Optional[str] = None,
        is_active: Optional[bool] = None,
        search: Optional[str] = None,
    ) -> List[PipelineTemplate]:
        """
        List pipeline templates with optional filters

        Args:
            category: Optional category filter
            is_active: Optional active status filter
            search: Optional search query (name, display_name, description)

        Returns:
            List of templates
        """
        query = self.db.query(PipelineTemplate)

        # Apply filters
        if category:
            query = query.filter(PipelineTemplate.category == category)

        if is_active is not None:
            query = query.filter(PipelineTemplate.is_active == is_active)

        if search:
            search_pattern = f"%{search}%"
            query = query.filter(
                or_(
                    PipelineTemplate.name.ilike(search_pattern),
                    PipelineTemplate.display_name.ilike(search_pattern),
                    PipelineTemplate.description.ilike(search_pattern)
                )
            )

        # Order by sort_order and display_name
        templates = query.order_by(
            PipelineTemplate.sort_order,
            PipelineTemplate.display_name
        ).all()

        return templates

    def get_categories(self) -> List[PipelineTemplateCategory]:
        """
        Get all pipeline categories with template counts

        Returns:
            List of categories with counts and template names
        """
        templates = self.db.query(PipelineTemplate).filter(
            PipelineTemplate.is_active == True
        ).all()

        # Group by category
        categories = defaultdict(lambda: {"count": 0, "templates": []})
        for template in templates:
            categories[template.category]["count"] += 1
            categories[template.category]["templates"].append(template.display_name)

        return [
            PipelineTemplateCategory(
                category=cat,
                count=data["count"],
                templates=data["templates"]
            )
            for cat, data in sorted(categories.items())
        ]

    # ============= CREATE OPERATIONS =============

    def create(self, template_data: PipelineTemplateCreate) -> PipelineTemplate:
        """
        Create a new custom pipeline template

        Args:
            template_data: Template creation data

        Returns:
            Created template

        Raises:
            HTTPException: If template name already exists
        """
        # Check if template name already exists
        existing = self.db.query(PipelineTemplate).filter(
            PipelineTemplate.name == template_data.name
        ).first()

        if existing:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Pipeline template with name '{template_data.name}' already exists"
            )

        # Create new template
        template = PipelineTemplate(
            **template_data.model_dump(),
            is_active=True,
            is_builtin=False  # Custom templates
        )

        self.db.add(template)
        self.db.commit()
        self.db.refresh(template)

        return template

    # ============= UPDATE OPERATIONS =============

    def update(
        self,
        template_id: UUID,
        update_data: PipelineTemplateUpdate
    ) -> PipelineTemplate:
        """
        Update pipeline template

        Args:
            template_id: Template UUID
            update_data: Update data

        Returns:
            Updated template

        Raises:
            HTTPException: If template not found

        Note:
            Built-in templates have restrictions on what can be updated
        """
        template = self.get_by_id_or_raise(template_id)

        # Update fields
        update_dict = update_data.model_dump(exclude_unset=True)

        # Prevent modifying critical fields for built-in templates
        if template.is_builtin:
            restricted_fields = ['name', 'script_name']
            for field in restricted_fields:
                if field in update_dict:
                    del update_dict[field]

        for field, value in update_dict.items():
            setattr(template, field, value)

        self.db.commit()
        self.db.refresh(template)

        return template

    def toggle_active(self, template_id: UUID) -> Tuple[PipelineTemplate, str]:
        """
        Toggle pipeline template active status

        Args:
            template_id: Template UUID

        Returns:
            Tuple of (updated template, status message)

        Raises:
            HTTPException: If template not found
        """
        template = self.get_by_id_or_raise(template_id)

        template.is_active = not template.is_active
        self.db.commit()

        status_text = "activated" if template.is_active else "deactivated"
        return template, f"Pipeline template {status_text} successfully"

    # ============= DELETE OPERATIONS =============

    def delete(self, template_id: UUID) -> None:
        """
        Delete pipeline template

        Args:
            template_id: Template UUID

        Raises:
            HTTPException: If template not found or is built-in
        """
        template = self.get_by_id_or_raise(template_id)

        if template.is_builtin:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Cannot delete built-in pipeline templates"
            )

        self.db.delete(template)
        self.db.commit()

    # ============= EXECUTION OPERATIONS =============

    def execute(
        self,
        user_id: UUID,
        execute_data: PipelineExecuteRequest
    ) -> PipelineTask:
        """
        Execute a pipeline with specified parameters

        Args:
            user_id: User UUID (for authorization)
            execute_data: Execution configuration

        Returns:
            Created pipeline task

        Raises:
            HTTPException: If template or project not found, or execution fails
        """
        # Verify template exists and is active
        template = self.db.query(PipelineTemplate).filter(
            PipelineTemplate.id == execute_data.template_id,
            PipelineTemplate.is_active == True
        ).first()

        if not template:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Pipeline template not found or inactive"
            )

        # Verify project belongs to user
        project = self.db.query(Project).filter(
            Project.id == execute_data.project_id,
            Project.user_id == user_id
        ).first()

        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Project not found"
            )

        # Merge default parameters with user-provided parameters
        merged_params = {**template.default_params, **execute_data.parameters}

        # Add template information to config
        config = {
            "template_id": str(execute_data.template_id),
            "template_name": template.name,
            "script_name": template.script_name,
            "script_path": template.script_path,
            "sample_ids": [str(sid) for sid in execute_data.sample_ids],
            "parameters": merged_params
        }

        # Create task
        task = PipelineTask(
            project_id=execute_data.project_id,
            task_name=execute_data.task_name,
            task_type=template.category,
            config=config,
            status="pending"
        )

        self.db.add(task)
        self.db.commit()
        self.db.refresh(task)

        # Submit to Celery for execution
        try:
            celery_task = run_ngs_pipeline.delay(
                task_id=str(task.id),
                pipeline_script=template.script_name,
                config=config
            )

            task.celery_task_id = celery_task.id
            task.status = "pending"
            task.started_at = datetime.utcnow()
            self.db.commit()
            self.db.refresh(task)

            return task

        except Exception as e:
            task.status = "failed"
            task.error_message = f"Failed to start pipeline: {str(e)}"
            self.db.commit()

            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Error starting pipeline: {str(e)}"
            )

    def batch_execute(
        self,
        user_id: UUID,
        execute_data: PipelineBatchExecuteRequest
    ) -> PipelineBatchExecuteResponse:
        """
        Execute a pipeline on multiple samples in batch

        Args:
            user_id: User UUID (for authorization)
            execute_data: Batch execution configuration

        Returns:
            Batch execution response with created tasks and failures

        Raises:
            HTTPException: If template, project, or samples not found
        """
        # Verify template exists and is active
        template = self.db.query(PipelineTemplate).filter(
            PipelineTemplate.id == execute_data.template_id,
            PipelineTemplate.is_active == True
        ).first()

        if not template:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Pipeline template not found or inactive"
            )

        # Verify project belongs to user
        project = self.db.query(Project).filter(
            Project.id == execute_data.project_id,
            Project.user_id == user_id
        ).first()

        if not project:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Project not found"
            )

        # Verify all samples belong to the project
        samples = self.db.query(Sample).filter(
            Sample.id.in_(execute_data.sample_ids),
            Sample.project_id == execute_data.project_id
        ).all()

        if len(samples) != len(execute_data.sample_ids):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Some samples not found or don't belong to the project"
            )

        # Merge default parameters with user-provided parameters
        merged_params = {**template.default_params, **execute_data.parameters}

        created_tasks = []
        failed_samples = []

        # Create and execute task for each sample
        for sample in samples:
            try:
                # Create task name with sample identifier
                task_name = f"{execute_data.task_name_prefix} - {sample.sample_id}"

                # Add template and sample information to config
                config = {
                    "template_id": str(execute_data.template_id),
                    "template_name": template.name,
                    "script_name": template.script_name,
                    "script_path": template.script_path,
                    "sample_ids": [str(sample.id)],
                    "sample_name": sample.sample_id,
                    "parameters": merged_params
                }

                # Create task
                task = PipelineTask(
                    project_id=execute_data.project_id,
                    task_name=task_name,
                    task_type=template.category,
                    config=config,
                    status="pending"
                )

                self.db.add(task)
                self.db.flush()  # Get task ID without committing

                # Submit to Celery for execution
                celery_task = run_ngs_pipeline.delay(
                    task_id=str(task.id),
                    pipeline_script=template.script_name,
                    config=config
                )

                task.celery_task_id = celery_task.id
                task.status = "pending"
                task.started_at = datetime.utcnow()

                created_tasks.append(task.id)

            except Exception as e:
                failed_samples.append({
                    "sample_id": str(sample.id),
                    "sample_name": sample.sample_id,
                    "error": str(e)
                })

        # Commit all successful tasks
        self.db.commit()

        return PipelineBatchExecuteResponse(
            total_tasks=len(created_tasks),
            created_tasks=created_tasks,
            failed_samples=failed_samples
        )

    # ============= AI RECOMMENDATION OPERATIONS =============

    def recommend_parameters(
        self,
        template_id: UUID,
        user_id: UUID,
        project_id: Optional[UUID] = None
    ) -> ParameterRecommendationResponse:
        """
        Get AI-powered parameter recommendations based on historical successful tasks

        Args:
            template_id: Template UUID
            user_id: User UUID (for authorization)
            project_id: Optional project ID for personalized recommendations

        Returns:
            Parameter recommendations with confidence score

        Raises:
            HTTPException: If template not found
        """
        # Verify template exists and is active
        template = self.db.query(PipelineTemplate).filter(
            PipelineTemplate.id == template_id,
            PipelineTemplate.is_active == True
        ).first()

        if not template:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Pipeline template not found or inactive"
            )

        # Query successful tasks for this template
        query = self.db.query(PipelineTask).filter(
            PipelineTask.config['template_id'].astext == str(template_id),
            PipelineTask.status == 'completed'
        )

        # If project_id provided, prioritize tasks from the same project
        if project_id:
            # Verify user has access to the project
            project = self.db.query(Project).filter(
                Project.id == project_id,
                Project.user_id == user_id
            ).first()

            if not project:
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND,
                    detail="Project not found"
                )

            # Get project-specific tasks first
            project_tasks = query.filter(PipelineTask.project_id == project_id).limit(20).all()
            if len(project_tasks) < 5:
                # If not enough project-specific tasks, include all user's tasks
                query = query.filter(PipelineTask.project_id.in_(
                    self.db.query(Project.id).filter(Project.user_id == user_id)
                ))
                all_tasks = query.limit(50).all()
            else:
                all_tasks = project_tasks
        else:
            # Get tasks from all user's projects
            query = query.filter(PipelineTask.project_id.in_(
                self.db.query(Project.id).filter(Project.user_id == user_id)
            ))
            all_tasks = query.limit(50).all()

        # If no historical tasks found, return template defaults
        if not all_tasks:
            return ParameterRecommendationResponse(
                recommended_params=template.default_params,
                confidence_score=0.5,
                based_on_tasks=0,
                explanation="No historical tasks found. Using template default parameters."
            )

        # Analyze parameters from successful tasks
        param_stats = defaultdict(lambda: defaultdict(int))
        total_tasks = len(all_tasks)

        for task in all_tasks:
            if 'parameters' in task.config:
                params = task.config['parameters']
                for key, value in params.items():
                    # Count frequency of each parameter value
                    param_stats[key][str(value)] += 1

        # Build recommended parameters
        recommended_params = {}
        for key, value_counts in param_stats.items():
            # Get most common value
            most_common_value = max(value_counts.items(), key=lambda x: x[1])
            value_str, count = most_common_value

            # Try to convert back to appropriate type
            try:
                if value_str.lower() in ['true', 'false']:
                    recommended_params[key] = value_str.lower() == 'true'
                elif '.' in value_str:
                    recommended_params[key] = float(value_str)
                elif value_str.isdigit():
                    recommended_params[key] = int(value_str)
                else:
                    recommended_params[key] = value_str
            except:
                recommended_params[key] = value_str

        # Merge with template defaults for missing parameters
        final_params = {**template.default_params, **recommended_params}

        # Calculate confidence score based on task count and parameter consistency
        if total_tasks >= 20:
            base_confidence = 0.9
        elif total_tasks >= 10:
            base_confidence = 0.8
        elif total_tasks >= 5:
            base_confidence = 0.7
        else:
            base_confidence = 0.6

        # Adjust confidence based on parameter consistency
        avg_consistency = sum(
            max(counts.values()) / sum(counts.values())
            for counts in param_stats.values()
        ) / len(param_stats) if param_stats else 0.5

        confidence_score = min(base_confidence * avg_consistency, 1.0)

        # Generate explanation
        if project_id and all_tasks[0].project_id == project_id:
            explanation = f"Recommendations based on {total_tasks} successful tasks from this project. "
        else:
            explanation = f"Recommendations based on {total_tasks} successful tasks from your projects. "

        explanation += f"These parameters were used in {int(avg_consistency * 100)}% of successful runs."

        return ParameterRecommendationResponse(
            recommended_params=final_params,
            confidence_score=round(confidence_score, 2),
            based_on_tasks=total_tasks,
            explanation=explanation
        )

    # ============= HELPER METHODS =============

    def get_active_templates_by_category(self, category: str) -> List[PipelineTemplate]:
        """
        Get active templates for a specific category

        Args:
            category: Pipeline category

        Returns:
            List of active templates in the category
        """
        return self.db.query(PipelineTemplate).filter(
            PipelineTemplate.category == category,
            PipelineTemplate.is_active == True
        ).order_by(
            PipelineTemplate.sort_order,
            PipelineTemplate.display_name
        ).all()

    def get_template_by_name(self, name: str) -> Optional[PipelineTemplate]:
        """
        Get template by name

        Args:
            name: Template name

        Returns:
            Template if found, None otherwise
        """
        return self.db.query(PipelineTemplate).filter(
            PipelineTemplate.name == name
        ).first()
