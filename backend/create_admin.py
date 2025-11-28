"""
Create default admin user

This script creates a default administrator account for first-time setup.
The admin password should be provided via environment variable ADMIN_PASSWORD.
"""
import asyncio
import logging
import os
from app.core.database import SessionLocal
from app.models.user import User
from app.core.security import get_password_hash
import uuid

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_admin_password() -> str:
    """
    Get admin password from environment variable.
    Falls back to a generated password if not set.
    """
    password = os.getenv("ADMIN_PASSWORD")
    if password:
        return password

    # Generate a random password if not provided
    import secrets
    generated = secrets.token_urlsafe(16)
    logger.warning("ADMIN_PASSWORD not set. Generated password: %s", generated)
    logger.warning("Please save this password and set ADMIN_PASSWORD environment variable for production.")
    return generated


async def create_admin_user():
    """
    Create default admin user
    """
    db = SessionLocal()

    try:
        # Check if admin already exists
        existing_admin = db.query(User).filter(User.username == "admin").first()

        if existing_admin:
            logger.info("[INFO] Admin user already exists")
            return

        admin_password = get_admin_password()

        # Create admin user
        admin_user = User(
            id=uuid.uuid4(),
            username="admin",
            email="admin@ngsmodule.com",
            password_hash=get_password_hash(admin_password),
            role="admin",
            is_active=True,
            storage_quota=1099511627776,  # 1TB for admin
            storage_used=0
        )

        db.add(admin_user)
        db.commit()

        logger.info("[OK] Default admin user created!")
        logger.info("   Username: admin")
        logger.info("[WARNING] Please change the password after first login!")

    except Exception as e:
        logger.error("[ERROR] Error creating admin user: %s", e)
        db.rollback()
    finally:
        db.close()


if __name__ == "__main__":
    asyncio.run(create_admin_user())
