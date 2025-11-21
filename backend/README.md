# NGSmodule Backend

FastAPI-based backend server for NGSmodule platform.

## Features

- ✅ RESTful API with FastAPI
- ✅ JWT Authentication
- ✅ PostgreSQL + SQLAlchemy ORM
- ✅ Alembic database migrations
- ✅ Celery for async tasks
- ✅ Redis for caching
- ✅ MinIO for file storage
- ✅ OpenAPI documentation

## Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Configure environment
cp .env.example .env
# Edit .env

# Run migrations
alembic upgrade head

# Start server
uvicorn app.main:app --reload
```

## API Documentation

Visit http://localhost:8000/api/v1/docs for interactive API documentation.

## Project Structure

```
app/
├── api/          # API routes
├── core/         # Configuration
├── models/       # Database models
├── schemas/      # Pydantic schemas
├── services/     # Business logic
└── main.py       # Application entry
```
