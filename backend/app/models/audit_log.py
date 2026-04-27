"""
Audit Log database model for tracking admin actions
"""
from sqlalchemy import Column, String, Text, DateTime, ForeignKey, JSON, Index
from sqlalchemy.orm import relationship
from app.core.types import UUID
import uuid

from app.core.database import Base
from app.core.datetime_utils import utc_now_naive


class AuditLog(Base):
    """
    Audit log model for tracking administrative actions

    Records all sensitive operations performed by administrators
    for compliance and security monitoring.
    """
    __tablename__ = "audit_logs"

    id = Column(UUID(), primary_key=True, default=uuid.uuid4)

    # Action information
    action = Column(String(100), nullable=False, index=True)
    timestamp = Column(DateTime, default=utc_now_naive, nullable=False, index=True)

    # Actor (admin who performed the action)
    admin_user_id = Column(UUID(), ForeignKey("users.id"), nullable=False, index=True)
    admin_username = Column(String(100), nullable=False)

    # Target (resource affected)
    target_user_id = Column(UUID(), ForeignKey("users.id"), nullable=True, index=True)
    target_username = Column(String(100), nullable=True)
    target_resource_type = Column(String(50), nullable=True)  # user, config, backup, etc.
    target_resource_id = Column(String(100), nullable=True)

    # Action details
    details = Column(JSON, default=dict)  # Additional context
    status = Column(String(20), default="success", index=True)  # success, failure

    # Request context
    ip_address = Column(String(45), nullable=True)  # IPv4 or IPv6
    user_agent = Column(Text, nullable=True)
    request_id = Column(String(100), nullable=True, index=True)

    # Relationships
    admin_user = relationship("User", foreign_keys=[admin_user_id])
    target_user = relationship("User", foreign_keys=[target_user_id])

    __table_args__ = (
        Index('ix_audit_logs_action_timestamp', 'action', 'timestamp'),
        Index('ix_audit_logs_admin_timestamp', 'admin_user_id', 'timestamp'),
    )

    def __repr__(self):
        return f"<AuditLog {self.action} by {self.admin_username} at {self.timestamp}>"

    def to_dict(self):
        """Convert audit log to dictionary"""
        return {
            "id": str(self.id),
            "action": self.action,
            "timestamp": self.timestamp.isoformat() if self.timestamp else None,
            "admin_user_id": str(self.admin_user_id),
            "admin_username": self.admin_username,
            "target_user_id": str(self.target_user_id) if self.target_user_id else None,
            "target_username": self.target_username,
            "target_resource_type": self.target_resource_type,
            "target_resource_id": self.target_resource_id,
            "details": self.details or {},
            "status": self.status,
            "ip_address": self.ip_address,
            "user_agent": self.user_agent,
            "request_id": self.request_id,
        }
