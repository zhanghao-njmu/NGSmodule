# 🧬 NGSmodule - 企业级生物信息学工作站

<div align="center">

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![React](https://img.shields.io/badge/react-18.2+-61dafb.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/fastapi-0.109+-009688.svg)](https://fastapi.tiangolo.com/)

**现代化的NGS数据处理平台 | 为研究人员和学生设计**

[快速开始](#-快速开始) •
[功能特性](#-功能特性) •
[技术栈](#-技术栈) •
[文档](#-文档)

</div>

---

## 📖 简介

NGSmodule 是一个企业级的生物信息学NGS数据处理平台，专为非编程背景的研究人员、学生和教师设计。

### 核心特性

- 🎯 **用户友好** - 现代化Web界面，无需编程经验
- ⚡ **高性能** - 异步任务处理，支持大规模数据分析
- 🔄 **实时监控** - WebSocket实时任务进度更新
- 🎨 **美观设计** - 基于Ant Design 5的现代UI
- 🚀 **易于部署** - Docker Compose一键启动

---

## 🚀 快速开始

### 一键启动

```bash
# 克隆项目
git clone <repository-url>
cd NGSmodule

# 运行快速启动脚本
./start.sh
```

### 访问应用

- 🌐 **Web界面**: http://localhost:3000
- 🔌 **API文档**: http://localhost:8000/api/v1/docs

默认登录: `admin` / `admin123`

详细部署说明: [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md)

---

## ✨ 功能特性

- ✅ **项目管理** - CRUD + 统计 + 归档
- ✅ **样本管理** - CSV批量导入
- ✅ **文件管理** - 50GB大文件 + MinIO存储
- ✅ **任务监控** - WebSocket实时进度

---

## 🏗️ 技术栈

**后端**: FastAPI + PostgreSQL + Redis + Celery + MinIO  
**前端**: React 18 + TypeScript + Ant Design + Zustand  
**部署**: Docker + Docker Compose

---

## 📊 开发进度

| 阶段 | 状态 | 完成度 |
|------|------|--------|
| Phase 1-3: 核心功能 | ✅ 完成 | 100% |
| Phase 4: 集成测试 | 🔄 进行中 | 80% |
| Phase 5-10: 高级功能 | 📋 计划中 | 0% |

**整体进度**: 30%

---

## 📚 文档

- [部署指南](DEPLOYMENT_GUIDE.md)
- [项目进度](PROJECT_PROGRESS_SUMMARY.md)
- [API文档](http://localhost:8000/api/v1/docs)

---

## 📄 许可证

MIT License

---

<div align="center">

**用 ❤️ 构建，为生物信息学研究服务**

</div>
