"""
Results API endpoints (Refactored with Service Layer)

This is the refactored version using ResultService for business logic.
Routers are now thin HTTP adapters that delegate to the service layer.
"""
from fastapi import APIRouter, Depends, Query
from sqlalchemy.orm import Session
from typing import Optional
from uuid import UUID

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.services.result_service import ResultService
from app.schemas.result import (
    ResultResponse,
    ResultListResponse,
    ResultVisualizationData,
)

router = APIRouter()


# ============= DEPENDENCY INJECTION =============

def get_result_service(db: Session = Depends(get_db)) -> ResultService:
    """Dependency to get ResultService instance"""
    return ResultService(db)


# ============= ENDPOINTS =============

@router.get("", response_model=ResultListResponse)
async def list_results(
    task_id: Optional[UUID] = Query(None, description="Filter by task ID"),
    result_type: Optional[str] = Query(None, description="Filter by result type"),
    skip: int = Query(0, ge=0),
    limit: int = Query(20, ge=1, le=100),
    current_user: User = Depends(get_current_user),
    service: ResultService = Depends(get_result_service)
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
    results, total = service.list_results(
        user_id=current_user.id,
        skip=skip,
        limit=limit,
        task_id=task_id,
        result_type=result_type
    )

    return ResultListResponse(
        results=results,
        total=total,
        skip=skip,
        limit=limit
    )


@router.get("/{result_id}", response_model=ResultResponse)
async def get_result(
    result_id: UUID,
    current_user: User = Depends(get_current_user),
    service: ResultService = Depends(get_result_service)
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
    return service.get_by_id_or_raise(result_id, current_user.id)


@router.get("/{result_id}/visualization", response_model=ResultVisualizationData)
async def get_visualization_data(
    result_id: UUID,
    current_user: User = Depends(get_current_user),
    service: ResultService = Depends(get_result_service)
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
    return service.get_visualization_data(result_id, current_user.id)


@router.get("/task/{task_id}/summary")
async def get_task_results_summary(
    task_id: UUID,
    current_user: User = Depends(get_current_user),
    service: ResultService = Depends(get_result_service)
):
    """
    Get a summary of all results for a specific task

    **Parameters:**
    - task_id: UUID of the task

    **Returns:**
    - Aggregated summary of all results associated with the task
    """
    return service.get_task_results_summary(task_id, current_user.id)
