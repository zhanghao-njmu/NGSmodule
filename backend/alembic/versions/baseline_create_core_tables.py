"""baseline: create the full Base.metadata schema

Revision ID: baseline_001
Revises:
Create Date: 2026-04-28 16:00:00.000000

This is the chain root. The chain previously had no migration creating
the core tables (users / projects / samples / ...), relying on
init_db.create_all running outside alembic — which made fresh
deployments fail with `no such table: users` once any later migration
tried to ALTER them.

Strategy: this baseline creates ALL current Base.metadata tables (14
total). The follow-up migrations (add_notifications_001 /
admin_enhanced_001 / data_downloads_001) are idempotent: they detect
that their target tables already exist and skip the body, so the chain
runs cleanly on fresh databases.

For legacy deployments where init_db + previous incremental ALTERs
already created the schema:
  alembic stamp data_downloads_001
to mark the DB as already at HEAD without re-running anything.

Future schema changes still need standard incremental ALTER migrations.
"""
from alembic import op


revision = 'baseline_001'
down_revision = None
branch_labels = None
depends_on = None


# Tables NOT created here (handled by later migrations):
#   notifications, notification_settings   → add_notifications_001
#   audit_logs, system_alerts, system_jobs, system_backups → admin_enhanced_001
#   download_jobs, vendor_credentials      → data_downloads_001


def upgrade():
    # Import inside the function so this file is safe to import before
    # the project is on sys.path (e.g. during alembic linting).
    from app.core.database import Base
    # Importing app.models triggers registration of every model on
    # Base.metadata; using create_all + checkfirst keeps the migration
    # idempotent and is portable across Postgres (prod) and SQLite (test).
    import app.models  # noqa: F401

    Base.metadata.create_all(bind=op.get_bind(), checkfirst=True)


def downgrade():
    from app.core.database import Base
    import app.models  # noqa: F401

    Base.metadata.drop_all(bind=op.get_bind(), checkfirst=True)
