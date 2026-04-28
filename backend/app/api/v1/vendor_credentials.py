"""
Vendor credential CRUD API.

Stores email/password for sequencing data-delivery vendors so users
don't have to retype on every download. Plaintext is never returned —
the frontend only sees masked previews.
"""
from typing import Optional
from uuid import UUID

from fastapi import APIRouter, Depends, Query, status
from sqlalchemy.orm import Session

from app.core.database import get_db
from app.core.deps import get_current_user
from app.models.user import User
from app.services.vendor_credential_service import VendorCredentialService
from app.schemas.vendor_credential import (
    VendorCredentialCreate,
    VendorCredentialResponse,
    VendorCredentialListResponse,
)
from app.schemas.common import MessageResponse

router = APIRouter()


def get_service(db: Session = Depends(get_db)) -> VendorCredentialService:
    return VendorCredentialService(db)


@router.post("", response_model=VendorCredentialResponse, status_code=status.HTTP_201_CREATED)
async def create_credential(
    payload: VendorCredentialCreate,
    user: User = Depends(get_current_user),
    service: VendorCredentialService = Depends(get_service),
):
    """Save (or update) a vendor credential. Same (vendor, label) pair is
    overwritten — equivalent to 'I changed my password'."""
    cred = service.create(user.id, payload)
    return VendorCredentialResponse(**service.to_response(cred))


@router.get("", response_model=VendorCredentialListResponse)
async def list_credentials(
    vendor: Optional[str] = Query(None, description="Filter by vendor id"),
    user: User = Depends(get_current_user),
    service: VendorCredentialService = Depends(get_service),
):
    """List the current user's saved credentials (masked)."""
    creds = service.list(user.id, vendor=vendor)
    items = [VendorCredentialResponse(**service.to_response(c)) for c in creds]
    return VendorCredentialListResponse(total=len(items), items=items)


@router.delete("/{cred_id}", response_model=MessageResponse)
async def delete_credential(
    cred_id: UUID,
    user: User = Depends(get_current_user),
    service: VendorCredentialService = Depends(get_service),
):
    """Remove a stored credential."""
    service.delete(cred_id, user.id)
    return MessageResponse(message="Credential deleted")
