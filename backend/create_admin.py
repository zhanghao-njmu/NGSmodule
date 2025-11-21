"""
Create default admin user

This script creates a default administrator account for first-time setup.
"""
import asyncio
from app.core.database import SessionLocal
from app.models.user import User
from app.core.security import get_password_hash
import uuid


async def create_admin_user():
    """
    Create default admin user
    """
    db = SessionLocal()

    try:
        # Check if admin already exists
        existing_admin = db.query(User).filter(User.username == "admin").first()

        if existing_admin:
            print("ℹ️  Admin user already exists")
            return

        # Create admin user
        admin_user = User(
            id=uuid.uuid4(),
            username="admin",
            email="admin@ngsmodule.com",
            password_hash=get_password_hash("admin123"),  # Default password
            role="admin",
            is_active=True,
            is_admin=True,
            storage_quota=1099511627776,  # 1TB for admin
            storage_used=0
        )

        db.add(admin_user)
        db.commit()

        print("✅ Default admin user created!")
        print("   Username: admin")
        print("   Password: admin123")
        print("   ⚠️  Please change the password after first login!")

    except Exception as e:
        print(f"❌ Error creating admin user: {e}")
        db.rollback()
    finally:
        db.close()


if __name__ == "__main__":
    asyncio.run(create_admin_user())
