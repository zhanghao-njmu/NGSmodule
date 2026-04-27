"""add admin enhanced tables (audit_logs, alerts, jobs, backups)

Revision ID: admin_enhanced_001
Revises: add_notifications_001
Create Date: 2026-04-27 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql


# revision identifiers, used by Alembic.
revision = 'admin_enhanced_001'
down_revision = 'add_notifications_001'
branch_labels = None
depends_on = None


def upgrade():
    # ====================================================================
    # audit_logs table
    # ====================================================================
    op.create_table(
        'audit_logs',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('action', sa.String(length=100), nullable=False),
        sa.Column('timestamp', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('admin_user_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('admin_username', sa.String(length=100), nullable=False),
        sa.Column('target_user_id', postgresql.UUID(as_uuid=True), nullable=True),
        sa.Column('target_username', sa.String(length=100), nullable=True),
        sa.Column('target_resource_type', sa.String(length=50), nullable=True),
        sa.Column('target_resource_id', sa.String(length=100), nullable=True),
        sa.Column('details', postgresql.JSON(astext_type=sa.Text()), nullable=True, server_default='{}'),
        sa.Column('status', sa.String(length=20), nullable=True, server_default='success'),
        sa.Column('ip_address', sa.String(length=45), nullable=True),
        sa.Column('user_agent', sa.Text(), nullable=True),
        sa.Column('request_id', sa.String(length=100), nullable=True),
        sa.ForeignKeyConstraint(['admin_user_id'], ['users.id']),
        sa.ForeignKeyConstraint(['target_user_id'], ['users.id']),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index('ix_audit_logs_action', 'audit_logs', ['action'])
    op.create_index('ix_audit_logs_timestamp', 'audit_logs', ['timestamp'])
    op.create_index('ix_audit_logs_admin_user_id', 'audit_logs', ['admin_user_id'])
    op.create_index('ix_audit_logs_target_user_id', 'audit_logs', ['target_user_id'])
    op.create_index('ix_audit_logs_status', 'audit_logs', ['status'])
    op.create_index('ix_audit_logs_request_id', 'audit_logs', ['request_id'])
    op.create_index('ix_audit_logs_action_timestamp', 'audit_logs', ['action', 'timestamp'])
    op.create_index('ix_audit_logs_admin_timestamp', 'audit_logs', ['admin_user_id', 'timestamp'])

    # ====================================================================
    # system_alerts table
    # ====================================================================
    op.create_table(
        'system_alerts',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('type', sa.String(length=20), nullable=False),
        sa.Column('severity', sa.String(length=20), nullable=False),
        sa.Column('title', sa.String(length=255), nullable=False),
        sa.Column('message', sa.Text(), nullable=False),
        sa.Column('source', sa.String(length=100), nullable=True),
        sa.Column('resolved', sa.Boolean(), nullable=False, server_default='false'),
        sa.Column('resolved_at', sa.DateTime(), nullable=True),
        sa.Column('resolved_by', postgresql.UUID(as_uuid=True), nullable=True),
        sa.Column('resolution_notes', sa.Text(), nullable=True),
        sa.Column('alert_metadata', postgresql.JSON(astext_type=sa.Text()), nullable=True, server_default='{}'),
        sa.Column('occurrence_count', sa.String(length=10), nullable=True, server_default='1'),
        sa.Column('timestamp', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('last_seen', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.ForeignKeyConstraint(['resolved_by'], ['users.id']),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index('ix_system_alerts_type', 'system_alerts', ['type'])
    op.create_index('ix_system_alerts_severity', 'system_alerts', ['severity'])
    op.create_index('ix_system_alerts_source', 'system_alerts', ['source'])
    op.create_index('ix_system_alerts_resolved', 'system_alerts', ['resolved'])
    op.create_index('ix_system_alerts_timestamp', 'system_alerts', ['timestamp'])
    op.create_index('ix_alerts_type_severity', 'system_alerts', ['type', 'severity'])
    op.create_index('ix_alerts_resolved_timestamp', 'system_alerts', ['resolved', 'timestamp'])

    # ====================================================================
    # system_jobs table
    # ====================================================================
    op.create_table(
        'system_jobs',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('type', sa.String(length=50), nullable=False),
        sa.Column('status', sa.String(length=20), nullable=False, server_default='pending'),
        sa.Column('user_id', postgresql.UUID(as_uuid=True), nullable=True),
        sa.Column('username', sa.String(length=100), nullable=True),
        sa.Column('celery_task_id', sa.String(length=255), nullable=True),
        sa.Column('progress', sa.Float(), nullable=True, server_default='0.0'),
        sa.Column('message', sa.Text(), nullable=True),
        sa.Column('parameters', postgresql.JSON(astext_type=sa.Text()), nullable=True, server_default='{}'),
        sa.Column('result', postgresql.JSON(astext_type=sa.Text()), nullable=True),
        sa.Column('error', sa.Text(), nullable=True),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('started_at', sa.DateTime(), nullable=True),
        sa.Column('completed_at', sa.DateTime(), nullable=True),
        sa.Column('retry_count', sa.String(length=10), nullable=True, server_default='0'),
        sa.Column('parent_job_id', postgresql.UUID(as_uuid=True), nullable=True),
        sa.ForeignKeyConstraint(['user_id'], ['users.id']),
        sa.ForeignKeyConstraint(['parent_job_id'], ['system_jobs.id']),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('celery_task_id')
    )
    op.create_index('ix_system_jobs_type', 'system_jobs', ['type'])
    op.create_index('ix_system_jobs_status', 'system_jobs', ['status'])
    op.create_index('ix_system_jobs_user_id', 'system_jobs', ['user_id'])
    op.create_index('ix_system_jobs_celery_task_id', 'system_jobs', ['celery_task_id'])
    op.create_index('ix_system_jobs_created_at', 'system_jobs', ['created_at'])
    op.create_index('ix_jobs_type_status', 'system_jobs', ['type', 'status'])
    op.create_index('ix_jobs_user_created', 'system_jobs', ['user_id', 'created_at'])

    # ====================================================================
    # system_backups table
    # ====================================================================
    op.create_table(
        'system_backups',
        sa.Column('id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('backup_type', sa.String(length=50), nullable=False),
        sa.Column('status', sa.String(length=20), nullable=False, server_default='pending'),
        sa.Column('file_path', sa.String(length=1024), nullable=True),
        sa.Column('file_name', sa.String(length=255), nullable=True),
        sa.Column('size', sa.BigInteger(), nullable=True, server_default='0'),
        sa.Column('compressed', sa.Boolean(), nullable=True, server_default='true'),
        sa.Column('checksum', sa.String(length=128), nullable=True),
        sa.Column('description', sa.Text(), nullable=True),
        sa.Column('backup_metadata', postgresql.JSON(astext_type=sa.Text()), nullable=True, server_default='{}'),
        sa.Column('error_message', sa.Text(), nullable=True),
        sa.Column('created_by', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.text('now()')),
        sa.Column('started_at', sa.DateTime(), nullable=True),
        sa.Column('completed_at', sa.DateTime(), nullable=True),
        sa.Column('expires_at', sa.DateTime(), nullable=True),
        sa.ForeignKeyConstraint(['created_by'], ['users.id']),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_index('ix_system_backups_backup_type', 'system_backups', ['backup_type'])
    op.create_index('ix_system_backups_status', 'system_backups', ['status'])
    op.create_index('ix_system_backups_created_at', 'system_backups', ['created_at'])
    op.create_index('ix_system_backups_expires_at', 'system_backups', ['expires_at'])
    op.create_index('ix_backups_type_status', 'system_backups', ['backup_type', 'status'])
    op.create_index('ix_backups_created_by_at', 'system_backups', ['created_by', 'created_at'])


def downgrade():
    # Drop system_backups
    op.drop_index('ix_backups_created_by_at', table_name='system_backups')
    op.drop_index('ix_backups_type_status', table_name='system_backups')
    op.drop_index('ix_system_backups_expires_at', table_name='system_backups')
    op.drop_index('ix_system_backups_created_at', table_name='system_backups')
    op.drop_index('ix_system_backups_status', table_name='system_backups')
    op.drop_index('ix_system_backups_backup_type', table_name='system_backups')
    op.drop_table('system_backups')

    # Drop system_jobs
    op.drop_index('ix_jobs_user_created', table_name='system_jobs')
    op.drop_index('ix_jobs_type_status', table_name='system_jobs')
    op.drop_index('ix_system_jobs_created_at', table_name='system_jobs')
    op.drop_index('ix_system_jobs_celery_task_id', table_name='system_jobs')
    op.drop_index('ix_system_jobs_user_id', table_name='system_jobs')
    op.drop_index('ix_system_jobs_status', table_name='system_jobs')
    op.drop_index('ix_system_jobs_type', table_name='system_jobs')
    op.drop_table('system_jobs')

    # Drop system_alerts
    op.drop_index('ix_alerts_resolved_timestamp', table_name='system_alerts')
    op.drop_index('ix_alerts_type_severity', table_name='system_alerts')
    op.drop_index('ix_system_alerts_timestamp', table_name='system_alerts')
    op.drop_index('ix_system_alerts_resolved', table_name='system_alerts')
    op.drop_index('ix_system_alerts_source', table_name='system_alerts')
    op.drop_index('ix_system_alerts_severity', table_name='system_alerts')
    op.drop_index('ix_system_alerts_type', table_name='system_alerts')
    op.drop_table('system_alerts')

    # Drop audit_logs
    op.drop_index('ix_audit_logs_admin_timestamp', table_name='audit_logs')
    op.drop_index('ix_audit_logs_action_timestamp', table_name='audit_logs')
    op.drop_index('ix_audit_logs_request_id', table_name='audit_logs')
    op.drop_index('ix_audit_logs_status', table_name='audit_logs')
    op.drop_index('ix_audit_logs_target_user_id', table_name='audit_logs')
    op.drop_index('ix_audit_logs_admin_user_id', table_name='audit_logs')
    op.drop_index('ix_audit_logs_timestamp', table_name='audit_logs')
    op.drop_index('ix_audit_logs_action', table_name='audit_logs')
    op.drop_table('audit_logs')
