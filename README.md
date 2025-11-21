# NGSmodule - Enterprise Bioinformatics Workstation

<div align="center">

![NGSmodule Logo](https://img.shields.io/badge/NGSmodule-v1.0.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.11-blue)
![React](https://img.shields.io/badge/react-18.2-blue)
![FastAPI](https://img.shields.io/badge/fastapi-0.109-green)

**A comprehensive, enterprise-grade platform for Next-Generation Sequencing data analysis**

[Features](#features) • [Quick Start](#quick-start) • [Documentation](#documentation) • [Architecture](#architecture)

</div>

---

## 🎯 Overview

NGSmodule is a modern, full-stack bioinformatics workstation designed for researchers and scientists to manage, analyze, and visualize NGS data through both **command-line interface** and an intuitive **web platform**.

### Key Highlights

- 🧬 **Comprehensive NGS Support**: RNA-seq, DNA-seq, scRNA-seq, ATAC-seq, ChIP-seq
- 🖥️ **Dual Interface**: Terminal CLI + Modern Web UI
- 🤖 **AI-Powered**: Intelligent parameter recommendations and result interpretation
- 🔐 **Enterprise Security**: JWT authentication, role-based access control
- 📊 **Rich Visualization**: Interactive charts with Plotly and ECharts
- 🐳 **Containerized**: Docker-based microservices architecture
- ⚡ **High Performance**: Async task processing with Celery
- 📱 **Responsive Design**: Works seamlessly on desktop, tablet, and mobile

---

## 🚀 Quick Start

### Prerequisites
- Docker & Docker Compose (recommended)

### Installation

```bash
# Clone repository
git clone https://github.com/zhanghao-njmu/NGSmodule.git
cd NGSmodule

# Copy environment files
cp backend/.env.example backend/.env
cp frontend/.env.example frontend/.env

# Start all services
docker-compose up -d

# Initialize database
docker-compose exec backend alembic upgrade head
```

### Access

- **Web UI**: http://localhost:3000
- **API Docs**: http://localhost:8000/api/v1/docs
- **Flower**: http://localhost:5555

## 📖 Documentation

- [Implementation Roadmap](IMPLEMENTATION_ROADMAP.md) - 10-phase plan
- [Progress Report](PROGRESS_REPORT.md) - Current status
- [Development Plan](DEVELOPMENT_PLAN.md) - Technical specs

## 📊 Status

**Overall Progress**: 15% (Phase 1 Complete)

| Phase | Status |
|-------|--------|
| Phase 1: Architecture | ✅ Complete |
| Phase 2: Backend Core | 🚧 Next |
| Phase 3: Frontend Core | 🚧 Next |

## 🤝 Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## 📝 License

MIT License - see [LICENSE](LICENSE) for details.

