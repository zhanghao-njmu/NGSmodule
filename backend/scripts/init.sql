-- NGSmodule Database Initialization Script
-- This script runs when the database is first created

-- Create extensions if needed
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";
CREATE EXTENSION IF NOT EXISTS "pg_trgm";  -- For text search optimization

-- Create custom types
DO $$ BEGIN
    CREATE TYPE project_status AS ENUM ('active', 'completed', 'archived');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

DO $$ BEGIN
    CREATE TYPE task_status AS ENUM ('pending', 'running', 'completed', 'failed', 'cancelled');
EXCEPTION
    WHEN duplicate_object THEN null;
END $$;

-- Enable row-level security (optional, for future use)
-- ALTER DATABASE ngsmodule SET row_security = on;

-- Set timezone
SET timezone = 'UTC';

-- Create indexes will be handled by Alembic migrations
-- This file is for database initialization only
