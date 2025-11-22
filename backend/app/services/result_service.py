"""
Result Service - Business logic for analysis result management
"""
from typing import List, Optional, Dict, Any, Tuple
from uuid import UUID
from sqlalchemy.orm import Session
from fastapi import HTTPException, status
import json
import random
import math
from pathlib import Path

from app.models.result import Result
from app.models.task import PipelineTask
from app.models.project import Project
from app.schemas.result import ResultVisualizationData


class ResultService:
    """Service class for result-related business logic"""

    def __init__(self, db: Session):
        """
        Initialize ResultService

        Args:
            db: Database session
        """
        self.db = db

    # ============= READ OPERATIONS =============

    def get_by_id(self, result_id: UUID, user_id: UUID) -> Optional[Result]:
        """
        Get result by ID for a specific user

        Args:
            result_id: Result UUID
            user_id: User UUID (for authorization)

        Returns:
            Result if found and belongs to user, None otherwise
        """
        result = self.db.query(Result).join(PipelineTask).join(Project).filter(
            Result.id == result_id,
            Project.user_id == user_id
        ).first()

        return result

    def get_by_id_or_raise(self, result_id: UUID, user_id: UUID) -> Result:
        """
        Get result by ID or raise 404 error

        Args:
            result_id: Result UUID
            user_id: User UUID (for authorization)

        Returns:
            Result

        Raises:
            HTTPException: If result not found or unauthorized
        """
        result = self.get_by_id(result_id, user_id)
        if not result:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Result {result_id} not found"
            )
        return result

    def list_results(
        self,
        user_id: UUID,
        skip: int = 0,
        limit: int = 20,
        task_id: Optional[UUID] = None,
        result_type: Optional[str] = None,
    ) -> Tuple[List[Result], int]:
        """
        List results for a user with pagination and filters

        Args:
            user_id: User UUID
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            task_id: Optional task filter
            result_type: Optional result type filter

        Returns:
            Tuple of (list of results, total count)
        """
        query = self.db.query(Result).join(PipelineTask).join(Project).filter(
            Project.user_id == user_id
        )

        # Apply filters
        if task_id:
            query = query.filter(Result.task_id == task_id)

        if result_type:
            query = query.filter(Result.result_type == result_type)

        # Get total count
        total = query.count()

        # Apply pagination
        results = query.order_by(Result.created_at.desc()).offset(skip).limit(limit).all()

        return results, total

    def get_task_results_summary(
        self,
        task_id: UUID,
        user_id: UUID
    ) -> Dict[str, Any]:
        """
        Get a summary of all results for a specific task

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)

        Returns:
            Aggregated summary of all results

        Raises:
            HTTPException: If task not found or unauthorized
        """
        # Verify task belongs to user
        task = self.db.query(PipelineTask).join(Project).filter(
            PipelineTask.id == task_id,
            Project.user_id == user_id
        ).first()

        if not task:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Task not found"
            )

        # Get results
        results = self.db.query(Result).filter(Result.task_id == task_id).all()

        if not results:
            return {
                "task_id": task_id,
                "total_results": 0,
                "result_types": [],
                "summary": "No results available yet",
            }

        # Aggregate summary
        result_types = list(set([r.result_type for r in results]))

        summary = {
            "task_id": task_id,
            "total_results": len(results),
            "result_types": result_types,
            "results_by_type": {},
        }

        for result_type in result_types:
            type_results = [r for r in results if r.result_type == result_type]
            summary["results_by_type"][result_type] = {
                "count": len(type_results),
                "latest": max(type_results, key=lambda x: x.created_at).created_at,
            }

        return summary

    # ============= VISUALIZATION OPERATIONS =============

    def get_visualization_data(
        self,
        result_id: UUID,
        user_id: UUID
    ) -> ResultVisualizationData:
        """
        Get visualization data for a specific result

        Args:
            result_id: Result UUID
            user_id: User UUID (for authorization)

        Returns:
            Structured data for visualization

        Raises:
            HTTPException: If result not found or unauthorized
        """
        result = self.get_by_id_or_raise(result_id, user_id)

        # Process result based on type
        if result.result_type == "qc_report":
            return self._generate_qc_visualization(result)
        elif result.result_type == "alignment":
            return self._generate_alignment_visualization(result)
        elif result.result_type == "quantification":
            return self._generate_quantification_visualization(result)
        elif result.result_type == "de_analysis":
            return self._generate_de_visualization(result)
        else:
            return ResultVisualizationData(
                type=result.result_type,
                message="Visualization not yet implemented for this result type",
                raw_metadata=result.metadata
            )

    # ============= HELPER METHODS =============

    def _generate_qc_visualization(self, result: Result) -> ResultVisualizationData:
        """
        Generate QC report visualization data

        Args:
            result: Result object

        Returns:
            QC visualization data
        """
        # In production, parse FastQC output from result.result_path
        # For now, use metadata or generate mock data
        metrics = {
            "total_reads": result.metadata.get("total_reads", 10000000),
            "quality_score": result.metadata.get("quality_score", 35.5),
            "gc_content": result.metadata.get("gc_content", 48.2),
            "duplication_rate": result.metadata.get("duplication_rate", 15.3),
        }

        charts = {
            "quality_distribution": {
                "x": list(range(1, 51)),  # Read positions
                "y": [30 + i * 0.2 for i in range(50)],  # Quality scores
                "type": "line",
            },
            "base_content": {
                "categories": ["A", "T", "G", "C", "N"],
                "values": [28.5, 28.3, 21.5, 21.5, 0.2],
                "type": "bar",
            },
            "gc_distribution": {
                "x": list(range(0, 101)),  # GC%
                "y": [abs(50 - i) * 10 for i in range(101)],  # Frequency
                "type": "area",
            },
        }

        status = "pass" if metrics["quality_score"] > 30 else "warning"

        return ResultVisualizationData(
            type="qc_report",
            metrics=metrics,
            charts=charts,
            status=status
        )

    def _generate_alignment_visualization(self, result: Result) -> ResultVisualizationData:
        """
        Generate alignment statistics visualization data

        Args:
            result: Result object

        Returns:
            Alignment visualization data
        """
        metrics = {
            "mapped_reads": result.metadata.get("mapped_reads", 8500000),
            "mapping_rate": result.metadata.get("mapping_rate", 85.0),
            "properly_paired": result.metadata.get("properly_paired", 7200000),
            "average_coverage": result.metadata.get("average_coverage", 42.5),
        }

        charts = {
            "mapping_summary": {
                "categories": ["Mapped", "Unmapped", "Duplicate"],
                "values": [85.0, 12.5, 2.5],
                "type": "pie",
            },
            "coverage_distribution": {
                "x": list(range(0, 101, 5)),  # Coverage bins
                "y": [100 - abs(50 - i) * 2 for i in range(0, 101, 5)],
                "type": "histogram",
            },
        }

        return ResultVisualizationData(
            type="alignment",
            metrics=metrics,
            charts=charts,
            status="pass"
        )

    def _generate_quantification_visualization(self, result: Result) -> ResultVisualizationData:
        """
        Generate gene quantification visualization data

        Args:
            result: Result object

        Returns:
            Quantification visualization data
        """
        # Generate mock gene expression data
        genes = [f"Gene_{i}" for i in range(100)]
        expression = [random.uniform(0, 1000) for _ in range(100)]

        metrics = {
            "total_genes": len(genes),
            "expressed_genes": sum(1 for e in expression if e > 10),
            "median_expression": sorted(expression)[len(expression) // 2],
        }

        charts = {
            "expression_distribution": {
                "x": genes[:20],  # Top 20 genes
                "y": sorted(expression, reverse=True)[:20],
                "type": "bar",
            },
            "density_plot": {
                "bins": 50,
                "data": expression,
                "type": "histogram",
            },
        }

        top_genes = [
            {"gene": genes[i], "expression": expression[i]}
            for i in sorted(range(len(expression)), key=lambda i: expression[i], reverse=True)[:10]
        ]

        return ResultVisualizationData(
            type="quantification",
            metrics=metrics,
            charts=charts,
            top_genes=top_genes
        )

    def _generate_de_visualization(self, result: Result) -> ResultVisualizationData:
        """
        Generate differential expression visualization data

        Args:
            result: Result object

        Returns:
            DE analysis visualization data
        """
        # Generate mock DE data
        genes = [f"Gene_{i}" for i in range(200)]
        log2fc = [random.uniform(-5, 5) for _ in range(200)]
        pvalues = [random.uniform(0, 1) for _ in range(200)]

        metrics = {
            "total_genes": len(genes),
            "up_regulated": sum(1 for fc in log2fc if fc > 1),
            "down_regulated": sum(1 for fc in log2fc if fc < -1),
            "significant": sum(1 for p in pvalues if p < 0.05),
        }

        charts = {
            "volcano_plot": {
                "x": log2fc,  # log2 fold change
                "y": [-math.log10(p) if p > 0 else 0 for p in pvalues],  # -log10 p-value
                "labels": genes,
                "type": "scatter",
            },
            "ma_plot": {
                "x": [random.uniform(0, 15) for _ in range(200)],  # A (average expression)
                "y": log2fc,  # M (log ratio)
                "type": "scatter",
            },
        }

        significant_genes = [
            {
                "gene": genes[i],
                "log2_fold_change": log2fc[i],
                "p_value": pvalues[i],
                "significant": pvalues[i] < 0.05 and abs(log2fc[i]) > 1,
            }
            for i in range(len(genes))
            if pvalues[i] < 0.05
        ][:20]  # Top 20 significant

        return ResultVisualizationData(
            type="de_analysis",
            metrics=metrics,
            charts=charts,
            significant_genes=significant_genes
        )

    def get_results_by_task(
        self,
        task_id: UUID,
        user_id: UUID,
        result_type: Optional[str] = None
    ) -> List[Result]:
        """
        Get all results for a specific task

        Args:
            task_id: Task UUID
            user_id: User UUID (for authorization)
            result_type: Optional result type filter

        Returns:
            List of results

        Raises:
            HTTPException: If task not found
        """
        # Verify task belongs to user
        task = self.db.query(PipelineTask).join(Project).filter(
            PipelineTask.id == task_id,
            Project.user_id == user_id
        ).first()

        if not task:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Task not found"
            )

        query = self.db.query(Result).filter(Result.task_id == task_id)

        if result_type:
            query = query.filter(Result.result_type == result_type)

        results = query.order_by(Result.created_at.desc()).all()

        return results
