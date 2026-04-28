"""
Common response schemas
"""

from typing import Any, Generic, List, Optional, TypeVar

from pydantic import BaseModel, Field


class MessageResponse(BaseModel):
    """Generic message response"""

    message: str
    detail: Optional[str] = None


class ErrorResponse(BaseModel):
    """Error response schema"""

    error: str
    detail: Optional[str] = None
    code: Optional[str] = None


class SuccessResponse(BaseModel):
    """Success response schema"""

    success: bool = True
    message: str
    data: Optional[Any] = None


# Generic type for pagination
T = TypeVar("T")


class PaginatedResponse(BaseModel, Generic[T]):
    """
    Unified paginated response schema
    统一的分页响应格式
    """

    total: int = Field(..., description="Total number of items")
    items: List[T] = Field(..., description="List of items")
    page: int = Field(..., description="Current page number (1-based)")
    page_size: int = Field(..., description="Number of items per page")
    has_more: bool = Field(False, description="Whether there are more pages")

    model_config = {"from_attributes": True}

    @classmethod
    def create(cls, items: List[T], total: int, page: int, page_size: int):
        """
        Create a paginated response

        Args:
            items: List of items for current page
            total: Total number of items
            page: Current page number (1-based)
            page_size: Number of items per page

        Returns:
            PaginatedResponse instance
        """
        return cls(total=total, items=items, page=page, page_size=page_size, has_more=(page * page_size) < total)


class PaginationParams(BaseModel):
    """
    Unified pagination parameters
    统一的分页参数
    """

    skip: int = Field(0, ge=0, description="Number of records to skip")
    limit: int = Field(20, ge=1, le=100, description="Number of records to return")

    @property
    def page(self) -> int:
        """Get current page number (1-based)"""
        return self.skip // self.limit + 1

    def apply(self, query):
        """
        Apply pagination to SQLAlchemy query

        Args:
            query: SQLAlchemy query object

        Returns:
            Query with offset and limit applied
        """
        return query.offset(self.skip).limit(self.limit)
