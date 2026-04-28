"""
User schemas for API request/response
"""
from pydantic import BaseModel, EmailStr, Field
from typing import Optional
from datetime import datetime
from uuid import UUID


class UserBase(BaseModel):
    """Base user schema"""
    username: str = Field(..., min_length=3, max_length=50)
    email: EmailStr
    full_name: Optional[str] = None
    organization: Optional[str] = None


class UserCreate(UserBase):
    """Schema for creating a new user"""
    password: str = Field(..., min_length=8)


class UserUpdate(BaseModel):
    """Schema for updating user"""
    full_name: Optional[str] = None
    organization: Optional[str] = None
    email: Optional[EmailStr] = None


class UserAdminUpdate(BaseModel):
    """Schema for admin updating user"""
    full_name: Optional[str] = None
    organization: Optional[str] = None
    email: Optional[EmailStr] = None
    role: Optional[str] = Field(None, pattern="^(user|admin)$")
    is_active: Optional[bool] = None
    storage_quota: Optional[int] = Field(None, ge=0)


class UserResponse(UserBase):
    """Schema for user response"""
    id: UUID
    role: str
    is_active: bool
    storage_quota: int
    storage_used: int
    created_at: datetime

    model_config = {"from_attributes": True}
class Token(BaseModel):
    """Schema for authentication token"""
    access_token: str
    token_type: str = "bearer"


class TokenPayload(BaseModel):
    """Schema for token payload"""
    sub: Optional[str] = None  # user_id


class UserStats(BaseModel):
    """Schema for user statistics"""
    user_id: UUID
    username: str
    total_projects: int
    total_samples: int
    total_tasks: int
    completed_tasks: int
    failed_tasks: int
    storage_used: int
    storage_quota: int
    storage_percent: float


class SystemStats(BaseModel):
    """Schema for system-wide statistics"""
    total_users: int
    active_users: int
    total_projects: int
    total_samples: int
    total_tasks: int
    running_tasks: int
    completed_tasks: int
    failed_tasks: int
    total_storage_used: int
    total_storage_quota: int
