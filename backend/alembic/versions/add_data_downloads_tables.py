"""add data_downloads tables (download_jobs, vendor_credentials)

Revision ID: data_downloads_001
Revises: admin_enhanced_001
Create Date: 2026-04-28 12:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


# revision identifiers, used by Alembic.
revision = 'data_downloads_001'
down_revision = 'admin_enhanced_001'
branch_labels = None
depends_on = None


def upgrade():
    # ====================================================================
    # download_jobs table — vendor data delivery tracking
    # ====================================================================
    op.create_table(
        'download_jobs',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('vendor', sa.String(length=32), nullable=False),
        sa.Column('source_path', sa.String(length=1024), nullable=False),
        sa.Column('dest_path', sa.String(length=1024), nullable=False),
        sa.Column('log_path', sa.String(length=1024), nullable=True),
        sa.Column('status', sa.String(length=16), nullable=False, server_default='pending'),
        sa.Column('progress_pct', sa.Float(), nullable=False, server_default='0.0'),
        sa.Column('bytes_downloaded', sa.BigInteger(), nullable=True, server_default='0'),
        sa.Column('file_size', sa.BigInteger(), nullable=True),
        sa.Column('error_message', sa.Text(), nullable=True),
        sa.Column('started_at', sa.DateTime(), nullable=True),
        sa.Column('finished_at', sa.DateTime(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        # Phase 5 auto-registration fields
        sa.Column('auto_register', sa.String(length=8), nullable=True),
        sa.Column('project_name_hint', sa.String(length=100), nullable=True),
        sa.Column('project_id', postgresql.UUID(as_uuid=True), nullable=True),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ondelete='CASCADE'),
        sa.ForeignKeyConstraint(['project_id'], ['projects.id'], ondelete='SET NULL'),
        sa.PrimaryKeyConstraint('id'),
    )
    op.create_index('ix_download_jobs_user_id', 'download_jobs', ['user_id'])
    op.create_index('ix_download_jobs_status', 'download_jobs', ['status'])
    op.create_index('ix_download_jobs_created_at', 'download_jobs', ['created_at'])

    # ====================================================================
    # vendor_credentials table — encrypted vendor login info
    # ====================================================================
    op.create_table(
        'vendor_credentials',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('vendor', sa.String(length=32), nullable=False),
        sa.Column('label', sa.String(length=64), nullable=False, server_default='default'),
        # Both encrypted_* columns store base64 fernet tokens (~150 chars typical, 512 covers paranoid).
        sa.Column('encrypted_email', sa.String(length=512), nullable=False),
        sa.Column('encrypted_password', sa.String(length=512), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('last_used_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['user_id'], ['users.id'], ondelete='CASCADE'),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('user_id', 'vendor', 'label', name='uq_vendor_cred_user_vendor_label'),
    )
    op.create_index('ix_vendor_credentials_user_id', 'vendor_credentials', ['user_id'])


def downgrade():
    op.drop_index('ix_vendor_credentials_user_id', table_name='vendor_credentials')
    op.drop_table('vendor_credentials')

    op.drop_index('ix_download_jobs_created_at', table_name='download_jobs')
    op.drop_index('ix_download_jobs_status', table_name='download_jobs')
    op.drop_index('ix_download_jobs_user_id', table_name='download_jobs')
    op.drop_table('download_jobs')
