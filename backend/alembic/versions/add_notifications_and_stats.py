"""add notifications and update user model

Revision ID: add_notifications_001
Revises:
Create Date: 2024-01-01 12:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

# revision identifiers, used by Alembic.
revision = 'add_notifications_001'
down_revision = 'baseline_001'  # chained after baseline (which creates users/projects/etc.)
branch_labels = None
depends_on = None


def upgrade():
    # Idempotent: baseline_001 creates the canonical Base.metadata schema
    # (which already includes last_login + the notifications/* tables),
    # so on a fresh deployment those objects exist. We guard each op so
    # the migration is also safe on a legacy DB that ran add_notifications
    # before baseline existed.
    bind = op.get_bind()
    inspector = sa.inspect(bind)

    if 'last_login' not in {c['name'] for c in inspector.get_columns('users')}:
        op.add_column('users', sa.Column('last_login', sa.DateTime(), nullable=True))

    if 'notifications' in inspector.get_table_names():
        return  # everything below already exists, baseline created it

    # Create notifications table
    op.create_table('notifications',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('type', sa.String(length=50), nullable=False),
        sa.Column('title', sa.String(length=255), nullable=False),
        sa.Column('message', sa.Text(), nullable=False),
        sa.Column('data', postgresql.JSON(astext_type=sa.Text()), nullable=True),
        sa.Column('read', sa.Boolean(), nullable=True, server_default='false'),
        sa.Column('action_url', sa.String(length=500), nullable=True),
        sa.Column('priority', sa.String(length=20), nullable=True, server_default='normal'),
        sa.Column('expires_at', sa.DateTime(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=True, server_default=sa.text('now()')),
        sa.Column('read_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index('ix_notifications_user_id', 'notifications', ['user_id'], unique=False)
    op.create_index('ix_notifications_type', 'notifications', ['type'], unique=False)
    op.create_index('ix_notifications_read', 'notifications', ['read'], unique=False)
    op.create_index('ix_notifications_created_at', 'notifications', ['created_at'], unique=False)

    # Create notification_settings table
    op.create_table('notification_settings',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('email_enabled', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('email_task_completed', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('email_task_failed', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('email_system_alerts', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('app_enabled', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('app_task_updates', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('app_project_updates', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('app_system_alerts', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('push_enabled', sa.Boolean(), nullable=True, server_default='false'),
        sa.Column('updated_at', sa.DateTime(), nullable=True, server_default=sa.text('now()')),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('user_id')
    )


def downgrade():
    # Drop notification_settings table
    op.drop_table('notification_settings')

    # Drop notifications table
    op.drop_index('ix_notifications_created_at', table_name='notifications')
    op.drop_index('ix_notifications_read', table_name='notifications')
    op.drop_index('ix_notifications_type', table_name='notifications')
    op.drop_index('ix_notifications_user_id', table_name='notifications')
    op.drop_table('notifications')

    # Remove last_login column from users table
    op.drop_column('users', 'last_login')
