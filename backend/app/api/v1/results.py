"""
Results API endpoints
"""
from typing import List, Optional
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from uuid import UUID

from app.core.database import get_db
from app.core.deps import get_current_user, verify_resource_ownership
from app.models.user import User
from app.models.result import Result
from app.models.task import PipelineTask
from app.schemas.result import (
    ResultResponse,
    ResultListResponse,
    ResultVisualizationData,
    QCMetrics,
    AlignmentStats,
)

router = APIRouter()


@router.get("", response_model=ResultListResponse)
async def list_results(
    task_id: Optional[UUID] = Query(None, description="Filter by task ID"),
    result_type: Optional[str] = Query(None, description="Filter by result type"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    """
    List results with optional filtering

    **Query Parameters:**
    - task_id: Filter results by specific task
    - result_type: Filter by result type (qc_report, alignment, quantification, de_analysis)
    - skip: Number of records to skip (pagination)
    - limit: Maximum number of records to return

    **Returns:**
    - List of results accessible to the current user
    """
    query = db.query(Result).join(PipelineTask)

    # Apply user filter (admin can see all, users see only their own)
    if current_user.role != "admin":
        query = query.filter(PipelineTask.user_id == current_user.id)

    # Apply filters
    if task_id:
        query = query.filter(Result.task_id == task_id)

    if result_type:
        query = query.filter(Result.result_type == result_type)

    # Get total count
    total = query.count()

    # Apply pagination
    results = query.order_by(Result.created_at.desc()).offset(skip).limit(limit).all()

    return {
        "results": results,
        "total": total,
        "skip": skip,
        "limit": limit,
    }


@router.get("/{result_id}", response_model=ResultResponse)
async def get_result(
    result_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    """
    Get a specific result by ID

    **Parameters:**
    - result_id: UUID of the result

    **Returns:**
    - Detailed result information

    **Raises:**
    - 404: Result not found
    - 403: User doesn't have permission to access this result
    """
    result = db.query(Result).filter(Result.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="Result not found")

    # Verify ownership
    verify_resource_ownership(db, PipelineTask, result.task_id, current_user)

    return result


@router.get("/{result_id}/visualization", response_model=ResultVisualizationData)
async def get_visualization_data(
    result_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    """
    Get visualization data for a specific result

    This endpoint processes the result files and returns data formatted
    for visualization in charts and graphs.

    **Parameters:**
    - result_id: UUID of the result

    **Returns:**
    - Structured data for visualization (charts, metrics, tables)

    **Supported Result Types:**
    - qc_report: Quality control metrics and charts
    - alignment: Alignment statistics and coverage plots
    - quantification: Expression levels and distributions
    - de_analysis: Differential expression results and volcano plots
    """
    result = db.query(Result).filter(Result.id == result_id).first()

    if not result:
        raise HTTPException(status_code=404, detail="Result not found")

    # Verify ownership
    verify_resource_ownership(db, PipelineTask, result.task_id, current_user)

    # Process result based on type
    visualization_data = _process_result_for_visualization(result)

    return visualization_data


@router.get("/task/{task_id}/summary")
async def get_task_results_summary(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
):
    """
    Get a summary of all results for a specific task

    **Parameters:**
    - task_id: UUID of the task

    **Returns:**
    - Aggregated summary of all results associated with the task
    """
    # Verify ownership
    verify_resource_ownership(db, PipelineTask, task_id, current_user)

    results = db.query(Result).filter(Result.task_id == task_id).all()

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


def _process_result_for_visualization(result: Result) -> dict:
    """
    Process result data for visualization

    This is a helper function that reads result files and formats them
    for frontend visualization.
    """
    import json
    from pathlib import Path

    # In a real implementation, you would:
    # 1. Read result files from result.result_path
    # 2. Parse the data (e.g., TSV, JSON, etc.)
    # 3. Format it for ECharts/Plotly

    # For now, return mock data based on result type
    if result.result_type == "qc_report":
        return _generate_qc_visualization(result)
    elif result.result_type == "alignment":
        return _generate_alignment_visualization(result)
    elif result.result_type == "quantification":
        return _generate_quantification_visualization(result)
    elif result.result_type == "de_analysis":
        return _generate_de_visualization(result)
    else:
        return {
            "type": result.result_type,
            "message": "Visualization not yet implemented for this result type",
            "raw_metadata": result.metadata,
        }


def _generate_qc_visualization(result: Result) -> dict:
    """Generate QC report visualization data"""
    # Mock data - in production, parse FastQC output
    return {
        "type": "qc_report",
        "metrics": {
            "total_reads": result.metadata.get("total_reads", 10000000),
            "quality_score": result.metadata.get("quality_score", 35.5),
            "gc_content": result.metadata.get("gc_content", 48.2),
            "duplication_rate": result.metadata.get("duplication_rate", 15.3),
        },
        "charts": {
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
        },
        "status": "pass" if result.metadata.get("quality_score", 35) > 30 else "warning",
    }


def _generate_alignment_visualization(result: Result) -> dict:
    """Generate alignment statistics visualization data"""
    return {
        "type": "alignment",
        "metrics": {
            "mapped_reads": result.metadata.get("mapped_reads", 8500000),
            "mapping_rate": result.metadata.get("mapping_rate", 85.0),
            "properly_paired": result.metadata.get("properly_paired", 7200000),
            "average_coverage": result.metadata.get("average_coverage", 42.5),
        },
        "charts": {
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
        },
        "status": "pass",
    }


def _generate_quantification_visualization(result: Result) -> dict:
    """Generate gene quantification visualization data"""
    import random

    # Generate mock gene expression data
    genes = [f"Gene_{i}" for i in range(100)]
    expression = [random.uniform(0, 1000) for _ in range(100)]

    return {
        "type": "quantification",
        "metrics": {
            "total_genes": len(genes),
            "expressed_genes": sum(1 for e in expression if e > 10),
            "median_expression": sorted(expression)[len(expression) // 2],
        },
        "charts": {
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
        },
        "top_genes": [
            {"gene": genes[i], "expression": expression[i]}
            for i in sorted(range(len(expression)), key=lambda i: expression[i], reverse=True)[:10]
        ],
    }


def _generate_de_visualization(result: Result) -> dict:
    """Generate differential expression visualization data"""
    import random
    import math

    # Generate mock DE data
    genes = [f"Gene_{i}" for i in range(200)]
    log2fc = [random.uniform(-5, 5) for _ in range(200)]
    pvalues = [random.uniform(0, 1) for _ in range(200)]

    return {
        "type": "de_analysis",
        "metrics": {
            "total_genes": len(genes),
            "up_regulated": sum(1 for fc in log2fc if fc > 1),
            "down_regulated": sum(1 for fc in log2fc if fc < -1),
            "significant": sum(1 for p in pvalues if p < 0.05),
        },
        "charts": {
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
        },
        "significant_genes": [
            {
                "gene": genes[i],
                "log2_fold_change": log2fc[i],
                "p_value": pvalues[i],
                "significant": pvalues[i] < 0.05 and abs(log2fc[i]) > 1,
            }
            for i in range(len(genes))
            if pvalues[i] < 0.05
        ][:20],  # Top 20 significant
    }
