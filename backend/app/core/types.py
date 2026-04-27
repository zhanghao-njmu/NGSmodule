"""
Cross-database type definitions.

Provides type aliases that gracefully degrade when running against SQLite
(used in unit tests) while preserving native PostgreSQL behavior in
production.
"""
from sqlalchemy import JSON, CHAR, TypeDecorator
from sqlalchemy.dialects.postgresql import UUID as PG_UUID, JSONB as PG_JSONB
import uuid as _uuid


class GUID(TypeDecorator):
    """Platform-independent UUID type.

    - PostgreSQL: native UUID
    - SQLite/others: CHAR(36) string

    The Python value is always a `uuid.UUID` (or string convertible to one).
    """

    impl = CHAR
    cache_ok = True

    def load_dialect_impl(self, dialect):
        if dialect.name == "postgresql":
            return dialect.type_descriptor(PG_UUID(as_uuid=True))
        return dialect.type_descriptor(CHAR(36))

    def process_bind_param(self, value, dialect):
        if value is None:
            return value
        if dialect.name == "postgresql":
            return value
        if not isinstance(value, _uuid.UUID):
            return str(_uuid.UUID(value))
        return str(value)

    def process_result_value(self, value, dialect):
        if value is None:
            return value
        if isinstance(value, _uuid.UUID):
            return value
        return _uuid.UUID(value)


class JSONType(TypeDecorator):
    """Platform-independent JSON type.

    - PostgreSQL: native JSONB
    - SQLite/others: generic JSON (stored as TEXT)
    """

    impl = JSON
    cache_ok = True

    def load_dialect_impl(self, dialect):
        if dialect.name == "postgresql":
            return dialect.type_descriptor(PG_JSONB())
        return dialect.type_descriptor(JSON())


# Convenience aliases — keep names matching the original PostgreSQL imports
# so model files only need to change the import path.
UUID = GUID
JSONB = JSONType


def UUIDColumn(*args, **kwargs):
    """Helper for UUID columns. Accepts the same `as_uuid` flag as
    `postgresql.UUID(as_uuid=True)` for backward compatibility but ignores it
    since GUID always returns a uuid.UUID."""
    kwargs.pop("as_uuid", None)
    return GUID(*args, **kwargs)
