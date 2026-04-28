"""
Vendor credential schemas.

Note: response schemas NEVER expose the raw email or password — only
metadata. To use a credential the client passes the credential id
when opening a vendor session; the backend decrypts internally.
"""
from datetime import datetime
from typing import Optional
from uuid import UUID

from pydantic import BaseModel, Field


class VendorCredentialCreate(BaseModel):
    vendor: str = Field(..., description="Vendor id, see /data-downloads/vendors")
    label: str = Field(default="default", max_length=64)
    email: str
    password: str = Field(..., min_length=1)


class VendorCredentialResponse(BaseModel):
    id: UUID
    vendor: str
    label: str
    email_preview: str  # masked for UI ('us***@example.com')
    created_at: datetime
    last_used_at: Optional[datetime] = None

    class Config:
        from_attributes = True


class VendorCredentialListResponse(BaseModel):
    total: int
    items: list[VendorCredentialResponse]


def mask_email(email: str) -> str:
    """Mask the local part of an email so the response can show a hint
    without leaking the full address."""
    if "@" not in email:
        return "***"
    local, domain = email.split("@", 1)
    if len(local) <= 2:
        masked = local[:1] + "*"
    else:
        masked = local[:2] + "*" * (len(local) - 2)
    return f"{masked}@{domain}"
