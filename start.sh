#!/bin/bash

# NGSmodule Quick Start Script
# This script helps you start the NGSmodule platform quickly

set -e

echo "🧬 NGSmodule Quick Start Script"
echo "================================"
echo ""

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check Docker
if ! command_exists docker; then
    echo "❌ Docker is not installed. Please install Docker first."
    echo "   Visit: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check Docker Compose
if ! command_exists docker-compose && ! docker compose version >/dev/null 2>&1; then
    echo "❌ Docker Compose is not installed. Please install Docker Compose first."
    echo "   Visit: https://docs.docker.com/compose/install/"
    exit 1
fi

echo "✅ Docker and Docker Compose are installed"
echo ""

# Create .env files if they don't exist
if [ ! -f backend/.env ]; then
    echo "📝 Creating backend/.env from example..."
    cp backend/.env.example backend/.env
fi

if [ ! -f frontend/.env ]; then
    echo "📝 Creating frontend/.env from example..."
    cp frontend/.env.example frontend/.env
fi

echo ""
echo "🚀 Starting NGSmodule services..."
echo "   This may take a few minutes on first run..."
echo ""

# Start services
docker-compose up -d

echo ""
echo "⏳ Waiting for services to be ready..."
sleep 10

# Check service health
echo ""
echo "🔍 Checking service health..."

# Check PostgreSQL
if docker-compose exec -T postgres pg_isready -U ngsmodule >/dev/null 2>&1; then
    echo "   ✅ PostgreSQL: Ready"
else
    echo "   ⚠️  PostgreSQL: Not ready yet"
fi

# Check Redis
if docker-compose exec -T redis redis-cli ping >/dev/null 2>&1; then
    echo "   ✅ Redis: Ready"
else
    echo "   ⚠️  Redis: Not ready yet"
fi

# Check Backend
if curl -s http://localhost:8000/health >/dev/null 2>&1; then
    echo "   ✅ Backend API: Ready"
else
    echo "   ⚠️  Backend API: Not ready yet (may need more time)"
fi

# Initialize database
echo ""
echo "🗄️  Initializing database..."
docker-compose exec -T backend python init_db.py

echo ""
echo "👤 Creating default admin user..."
docker-compose exec -T backend python create_admin.py

echo ""
echo "✨ NGSmodule is ready!"
echo ""
echo "📍 Access the application:"
echo "   🌐 Frontend:     http://localhost:3000"
echo "   🔌 Backend API:  http://localhost:8000"
echo "   📚 API Docs:     http://localhost:8000/api/v1/docs"
echo "   🌺 Flower:       http://localhost:5555"
echo "   💾 MinIO:        http://localhost:9001"
echo ""
echo "👤 Default credentials:"
echo "   Username: admin"
echo "   Password: admin123"
echo "   ⚠️  Please change the password after first login!"
echo ""
echo "📝 View logs:"
echo "   docker-compose logs -f"
echo ""
echo "🛑 Stop services:"
echo "   docker-compose down"
echo ""
