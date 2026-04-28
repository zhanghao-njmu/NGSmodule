"""
Vendor credential model — encrypted email/password for sequencing
data-delivery vendors (联川 / 诺禾致源 / ...).

Stored ciphertexts are fernet-encrypted with a key derived from
`settings.SECRET_KEY`. Rotating SECRET_KEY invalidates stored
credentials (intentional — forces re-entry).
"""

import uuid

from sqlalchemy import Column, DateTime, ForeignKey, String, UniqueConstraint
from sqlalchemy.orm import relationship

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive
from app.core.types import UUID


class VendorCredential(Base):
    __tablename__ = "vendor_credentials"
    __table_args__ = (UniqueConstraint("user_id", "vendor", "label", name="uq_vendor_cred_user_vendor_label"),)

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(), ForeignKey("users.id", ondelete="CASCADE"), nullable=False)
    vendor = Column(String(32), nullable=False)
    # Optional label distinguishes multiple credentials for the same vendor
    # (e.g. 'lab pi account' vs 'shared bioinfo account').
    label = Column(String(64), nullable=False, default="default")

    # Both fields hold base64 fernet tokens — never plaintext.
    encrypted_email = Column(String(512), nullable=False)
    encrypted_password = Column(String(512), nullable=False)

    created_at = Column(DateTime, default=utc_now_naive, nullable=False)
    last_used_at = Column(DateTime)

    user = relationship("User")

    def __repr__(self):
        return f"<VendorCredential vendor={self.vendor} label={self.label}>"
