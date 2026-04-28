"""
Vendor credential service — CRUD + decrypt-on-demand.
"""
from typing import List, Optional, Tuple
from uuid import UUID

from fastapi import HTTPException, status
from sqlalchemy.orm import Session

from app.core.datetime_utils import utc_now_naive
from app.core.security import encrypt_secret, decrypt_secret
from app.models.vendor_credential import VendorCredential
from app.schemas.vendor_credential import VendorCredentialCreate, mask_email


class VendorCredentialService:
    def __init__(self, db: Session):
        self.db = db

    def create(self, user_id: UUID, payload: VendorCredentialCreate) -> VendorCredential:
        # Re-creating an existing (vendor, label) overwrites — UX matches
        # 'I'm updating my password'.
        existing = (
            self.db.query(VendorCredential)
            .filter(
                VendorCredential.user_id == user_id,
                VendorCredential.vendor == payload.vendor,
                VendorCredential.label == payload.label,
            )
            .first()
        )
        if existing:
            existing.encrypted_email = encrypt_secret(payload.email)
            existing.encrypted_password = encrypt_secret(payload.password)
            self.db.commit()
            self.db.refresh(existing)
            return existing

        cred = VendorCredential(
            user_id=user_id,
            vendor=payload.vendor,
            label=payload.label,
            encrypted_email=encrypt_secret(payload.email),
            encrypted_password=encrypt_secret(payload.password),
        )
        self.db.add(cred)
        self.db.commit()
        self.db.refresh(cred)
        return cred

    def list(self, user_id: UUID, vendor: Optional[str] = None) -> List[VendorCredential]:
        q = self.db.query(VendorCredential).filter(VendorCredential.user_id == user_id)
        if vendor:
            q = q.filter(VendorCredential.vendor == vendor)
        return q.order_by(VendorCredential.created_at.desc()).all()

    def get(self, cred_id: UUID, user_id: UUID) -> VendorCredential:
        cred = (
            self.db.query(VendorCredential)
            .filter(VendorCredential.id == cred_id, VendorCredential.user_id == user_id)
            .first()
        )
        if not cred:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Vendor credential {cred_id} not found",
            )
        return cred

    def delete(self, cred_id: UUID, user_id: UUID) -> None:
        cred = self.get(cred_id, user_id)
        self.db.delete(cred)
        self.db.commit()

    def reveal(self, cred_id: UUID, user_id: UUID) -> Tuple[str, str]:
        """Decrypt and mark as used. Internal use only — never exposed via API."""
        cred = self.get(cred_id, user_id)
        try:
            email = decrypt_secret(cred.encrypted_email)
            password = decrypt_secret(cred.encrypted_password)
        except Exception as exc:
            raise HTTPException(
                status_code=status.HTTP_410_GONE,
                detail=(
                    "Stored credential could not be decrypted "
                    "(SECRET_KEY may have rotated); please re-enter."
                ),
            ) from exc
        cred.last_used_at = utc_now_naive()
        self.db.commit()
        return email, password

    def to_response(self, cred: VendorCredential) -> dict:
        """Build the response dict (with masked email) for API."""
        # Best-effort decrypt for masking; if it fails we still return the
        # row so user can delete/recreate.
        try:
            email = decrypt_secret(cred.encrypted_email)
        except Exception:
            email = "(decryption-failed)"
        return {
            "id": cred.id,
            "vendor": cred.vendor,
            "label": cred.label,
            "email_preview": mask_email(email),
            "created_at": cred.created_at,
            "last_used_at": cred.last_used_at,
        }
