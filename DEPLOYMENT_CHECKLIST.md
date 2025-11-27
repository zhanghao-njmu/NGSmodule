# NGSmodule 生产部署检查清单

**项目**: NGSmodule 生物信息学工作站
**版本**: 1.0.0
**更新日期**: 2024-01-01

---

## 📋 使用说明

本清单用于确保 NGSmodule 项目在部署到生产环境前已完成所有必要的准备工作。

**使用方法**:
1. 逐项检查每个条目
2. 完成后在对应的 `[ ]` 中标记 `[x]`
3. 确保所有必需项（标记为 ⚠️）都已完成
4. 可选项（标记为 💡）根据实际需求决定是否完成

---

## 🎯 Phase 27 完成验证

### 代码质量

- [x] TypeScript编译无错误 (`npm run build`)
- [x] ESLint检查通过
- [x] 代码冗余已消除
- [x] 废弃类型已删除
- [x] 工具函数使用统一

### 性能优化

- [x] 路由懒加载已实现
- [x] 代码分割完成（15个chunks）
- [x] Bundle大小优化（-55%）
- [x] Gzip压缩配置
- [x] 缓存策略优化

### 系统一致性

- [x] API服务层使用统一模式
- [x] 状态管理模式一致
- [x] TypeScript类型系统完整
- [x] UI组件库统一（Ant Design 5）

### 生产配置

- [x] Nginx配置增强完成
- [x] 反向代理配置创建
- [x] Docker生产编排完成
- [x] 环境变量模板创建
- [x] 自动化部署脚本创建
- [x] 生产部署文档完成

---

## 🔧 环境准备

### 服务器环境 ⚠️

- [ ] 服务器已准备（满足最低配置要求）
  - CPU: 4核心+
  - 内存: 8GB+
  - 存储: 100GB+ SSD
  - 网络: 100Mbps+

- [ ] 操作系统已更新
  ```bash
  sudo apt update && sudo apt upgrade -y
  ```

- [ ] Docker已安装并运行
  ```bash
  docker --version  # 应≥20.10.0
  docker ps  # 验证Docker运行正常
  ```

- [ ] Docker Compose已安装
  ```bash
  docker-compose --version  # 应≥2.0.0
  ```

- [ ] Git已安装
  ```bash
  git --version  # 应≥2.25.0
  ```

- [ ] OpenSSL已安装
  ```bash
  openssl version  # 应≥1.1.1
  ```

### 网络配置 ⚠️

- [ ] 服务器有公网IP地址
- [ ] 域名已注册（如需使用域名）
- [ ] DNS记录已配置
  - A记录: yourdomain.com → 服务器IP
  - A记录: www.yourdomain.com → 服务器IP
  - CNAME (可选): api.yourdomain.com → yourdomain.com

- [ ] 防火墙规则已配置
  ```bash
  # 必需端口
  sudo ufw allow 22/tcp   # SSH
  sudo ufw allow 80/tcp   # HTTP
  sudo ufw allow 443/tcp  # HTTPS
  sudo ufw enable
  ```

- [ ] 不必要的端口已关闭
  - 5432 (PostgreSQL) - 应仅内部访问
  - 6379 (Redis) - 应仅内部访问
  - 9000/9001 (MinIO) - 根据需求决定

---

## 🔐 安全配置

### SSL/TLS证书 ⚠️

**选项A: Let's Encrypt（推荐生产环境）**
- [ ] Certbot已安装
  ```bash
  sudo apt install certbot
  ```

- [ ] SSL证书已获取
  ```bash
  sudo certbot certonly --standalone -d yourdomain.com -d www.yourdomain.com
  ```

- [ ] 证书已复制到项目目录
  ```bash
  sudo cp /etc/letsencrypt/live/yourdomain.com/fullchain.pem ./ssl/cert.pem
  sudo cp /etc/letsencrypt/live/yourdomain.com/privkey.pem ./ssl/key.pem
  sudo chown $USER:$USER ./ssl/*.pem
  ```

- [ ] 证书自动续期已配置
  ```bash
  sudo crontab -e
  # 添加: 0 3 * * * certbot renew --quiet
  ```

**选项B: 自签名证书（仅测试环境）**
- [ ] 自签名证书已生成
  ```bash
  openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
    -keyout ssl/key.pem -out ssl/cert.pem
  ```

### 密钥和密码 ⚠️

- [ ] 所有默认密码已更改
- [ ] 强随机密钥已生成
  ```bash
  # SECRET_KEY
  openssl rand -hex 32

  # JWT_SECRET_KEY
  openssl rand -hex 32

  # POSTGRES_PASSWORD
  openssl rand -base64 32

  # REDIS_PASSWORD
  openssl rand -base64 32

  # MINIO_ROOT_PASSWORD
  openssl rand -base64 32

  # FLOWER_PASSWORD
  openssl rand -base64 16
  ```

- [ ] `.env` 文件权限正确
  ```bash
  chmod 600 .env
  ```

- [ ] `.env` 文件已添加到 `.gitignore`
- [ ] 敏感信息未硬编码在代码中

---

## ⚙️ 配置文件

### 环境变量 ⚠️

- [ ] `.env` 文件已从模板创建
  ```bash
  cp .env.production .env
  ```

- [ ] 所有必需变量已配置:
  - [ ] `SECRET_KEY` (64位hex字符串)
  - [ ] `JWT_SECRET_KEY` (64位hex字符串)
  - [ ] `POSTGRES_PASSWORD` (强密码)
  - [ ] `REDIS_PASSWORD` (强密码)
  - [ ] `MINIO_ROOT_USER` (用户名)
  - [ ] `MINIO_ROOT_PASSWORD` (强密码)
  - [ ] `FLOWER_PASSWORD` (强密码)

- [ ] 域名配置已更新:
  - [ ] `VITE_API_URL` = `https://yourdomain.com/api/v1` 或 `https://api.yourdomain.com/api/v1`
  - [ ] `VITE_WS_URL` = `wss://yourdomain.com/ws` 或 `wss://api.yourdomain.com/ws`
  - [ ] `BACKEND_CORS_ORIGINS` = `https://yourdomain.com,https://www.yourdomain.com`

- [ ] 性能配置已调整（根据服务器规格）:
  - [ ] `UVICORN_WORKERS` (推荐: CPU核心数)
  - [ ] `CELERY_WORKERS` (推荐: 2-4)
  - [ ] `CELERY_CONCURRENCY` (推荐: 4-8)

### 💡 可选配置

- [ ] 邮件服务已配置（如需通知功能）:
  - [ ] `SMTP_HOST`
  - [ ] `SMTP_USER`
  - [ ] `SMTP_PASSWORD`

- [ ] Sentry已配置（如需错误追踪）:
  - [ ] `SENTRY_DSN`

- [ ] 备份策略已配置:
  - [ ] `BACKUP_ENABLED=true`
  - [ ] `BACKUP_SCHEDULE` (cron表达式)
  - [ ] `BACKUP_RETENTION_DAYS`

### Docker配置

- [ ] `docker-compose.prod.yml` 已审查
- [ ] 服务资源限制已配置（如需要）
- [ ] 持久化卷路径已确认
- [ ] 网络配置已确认

### Nginx配置

- [ ] `nginx-proxy.conf` 中的域名已更新
  ```nginx
  server_name yourdomain.com;
  ```

- [ ] SSL配置已启用（如使用HTTPS）
  ```nginx
  listen 443 ssl http2;
  ssl_certificate /etc/nginx/ssl/cert.pem;
  ssl_certificate_key /etc/nginx/ssl/key.pem;
  ```

- [ ] 上传文件大小限制已确认
  ```nginx
  client_max_body_size 500M;
  ```

---

## 🚀 部署执行

### 代码准备

- [ ] 代码已克隆到服务器
  ```bash
  git clone <repository-url>
  cd NGSmodule
  ```

- [ ] 分支已切换到生产分支
  ```bash
  git checkout main  # 或你的生产分支
  ```

- [ ] 代码是最新版本
  ```bash
  git pull origin main
  ```

### 目录和权限

- [ ] 必需目录已创建
  ```bash
  mkdir -p backups/postgres ssl logs
  ```

- [ ] 目录权限已设置
  ```bash
  chmod 700 backups ssl logs
  chmod 600 .env
  ```

- [ ] SSL证书权限正确
  ```bash
  chmod 600 ssl/*.pem
  ```

### 初始化部署 ⚠️

**选项A: 使用自动化脚本（推荐）**

- [ ] 部署脚本可执行
  ```bash
  chmod +x deploy-production.sh
  ```

- [ ] 运行初始化设置
  ```bash
  ./deploy-production.sh setup
  ```

- [ ] 编辑并验证 `.env` 配置
  ```bash
  nano .env
  ```

- [ ] 构建Docker镜像
  ```bash
  ./deploy-production.sh build
  ```

- [ ] 启动所有服务
  ```bash
  ./deploy-production.sh start
  ```

- [ ] 检查服务状态
  ```bash
  ./deploy-production.sh status
  ```

**选项B: 手动部署**

- [ ] 构建Docker镜像
  ```bash
  docker-compose -f docker-compose.prod.yml build --no-cache
  ```

- [ ] 启动所有服务
  ```bash
  docker-compose -f docker-compose.prod.yml up -d
  ```

- [ ] 查看服务状态
  ```bash
  docker-compose -f docker-compose.prod.yml ps
  ```

- [ ] 查看日志
  ```bash
  docker-compose -f docker-compose.prod.yml logs -f
  ```

---

## ✅ 部署验证

### 服务健康检查 ⚠️

- [ ] 所有容器正在运行
  ```bash
  docker-compose -f docker-compose.prod.yml ps
  # 所有服务应显示 "Up" 状态
  ```

- [ ] PostgreSQL健康检查通过
  ```bash
  docker-compose -f docker-compose.prod.yml exec postgres pg_isready
  # 应显示: accepting connections
  ```

- [ ] Redis健康检查通过
  ```bash
  docker-compose -f docker-compose.prod.yml exec redis redis-cli ping
  # 应返回: PONG
  ```

- [ ] 后端API健康检查通过
  ```bash
  curl http://localhost:8000/health
  # 或
  curl https://yourdomain.com/api/health
  # 应返回: {"status":"healthy"}
  ```

- [ ] 前端可访问
  ```bash
  curl http://localhost/
  # 或
  curl https://yourdomain.com/
  # 应返回HTML内容
  ```

- [ ] Nginx状态正常
  ```bash
  curl http://localhost:8080/nginx-status
  # 应显示nginx状态信息
  ```

### 功能测试 ⚠️

- [ ] 可以访问首页
  - 浏览器打开: `https://yourdomain.com`

- [ ] 登录功能正常
  - 测试登录页面
  - 测试注册页面（如果启用）

- [ ] API端点可访问
  ```bash
  curl https://yourdomain.com/api/v1/health
  ```

- [ ] WebSocket连接正常（如果使用）
  - 测试实时功能
  - 检查浏览器控制台无WebSocket错误

- [ ] 静态资源加载正常
  - 检查浏览器Network标签
  - CSS、JS、图片等资源正常加载
  - 无404错误

### 性能验证

- [ ] 页面加载时间acceptable
  - 首屏加载 < 3秒 (4G网络)
  - 页面切换流畅

- [ ] 资源大小合理
  - 主bundle < 100KB
  - Gzipped < 40KB

- [ ] 缓存正常工作
  - 静态资源有缓存headers
  - 重复访问使用缓存

### 安全验证 ⚠️

- [ ] HTTPS正常工作（如果配置）
  ```bash
  curl -I https://yourdomain.com
  # 检查返回200且使用HTTPS
  ```

- [ ] 安全headers存在
  ```bash
  curl -I https://yourdomain.com | grep -i "x-frame-options\|x-content-type-options\|x-xss-protection"
  ```

- [ ] 敏感端口不对外开放
  ```bash
  # 从外部测试（不是服务器本地）
  telnet yourdomain.com 5432  # 应连接失败
  telnet yourdomain.com 6379  # 应连接失败
  ```

- [ ] 敏感文件无法访问
  ```bash
  curl https://yourdomain.com/.env  # 应返回403或404
  curl https://yourdomain.com/.git/config  # 应返回403或404
  ```

---

## 📊 监控和日志

### 日志配置

- [ ] 日志目录已创建
  ```bash
  mkdir -p logs
  ```

- [ ] 日志轮转正常工作
  - 检查Docker日志配置
  - 确认日志文件大小限制

- [ ] 可以查看实时日志
  ```bash
  ./deploy-production.sh logs
  # 或
  docker-compose -f docker-compose.prod.yml logs -f
  ```

### 监控设置 💡

- [ ] Flower (Celery监控) 可访问
  - 访问: `http://yourdomain.com:5555`
  - 使用配置的用户名和密码登录

- [ ] Nginx状态监控可访问
  ```bash
  curl http://localhost:8080/nginx-status
  ```

- [ ] 系统资源监控已设置（可选）
  ```bash
  docker stats
  ```

- [ ] 错误追踪已配置（如使用Sentry）
  - Sentry DSN已配置
  - 测试错误上报

---

## 💾 备份配置

### 备份策略 ⚠️

- [ ] 备份目录已创建
  ```bash
  mkdir -p backups/postgres
  chmod 700 backups
  ```

- [ ] 可以手动备份数据库
  ```bash
  ./deploy-production.sh backup
  ```

- [ ] 备份文件生成成功
  ```bash
  ls -lh backups/postgres/
  ```

- [ ] 自动备份已配置
  ```bash
  # 配置cron任务
  crontab -e
  # 添加: 0 2 * * * cd /path/to/NGSmodule && ./deploy-production.sh backup
  ```

- [ ] 备份恢复已测试
  ```bash
  # 在测试环境测试恢复
  ./deploy-production.sh restore backups/postgres/backup_YYYYMMDD_HHMMSS.sql.gz
  ```

### 备份内容清单

- [ ] 数据库备份 (PostgreSQL)
- [ ] 配置文件备份 (.env)
- [ ] SSL证书备份 (ssl/*)
- [ ] 上传文件备份（根据需要）
- [ ] 重要日志备份（根据需要）

---

## 🔄 运维准备

### 文档 ⚠️

- [ ] 生产部署文档已阅读
  - `PRODUCTION_DEPLOYMENT.md`

- [ ] Phase 27报告已阅读
  - `PHASE_27_COMPLETE_REPORT.md`

- [ ] 故障排查指南已熟悉
  - 查阅 `PRODUCTION_DEPLOYMENT.md` 故障排查章节

### 运维流程

- [ ] 知道如何查看日志
  ```bash
  ./deploy-production.sh logs
  ```

- [ ] 知道如何重启服务
  ```bash
  ./deploy-production.sh restart
  ```

- [ ] 知道如何更新部署
  ```bash
  ./deploy-production.sh update
  ```

- [ ] 知道如何备份和恢复
  ```bash
  ./deploy-production.sh backup
  ./deploy-production.sh restore <backup-file>
  ```

### 联系人信息

- [ ] 技术负责人联系方式已记录
- [ ] 紧急联系流程已建立
- [ ] 相关文档位置已记录

---

## 📝 最终检查

### 部署前最终确认 ⚠️

- [ ] 所有必需项 (标记⚠️) 已完成
- [ ] 配置文件已double-check
- [ ] 密钥和密码已更新并记录（安全存储）
- [ ] DNS已生效（如使用域名）
- [ ] SSL证书有效（如使用HTTPS）
- [ ] 备份策略已测试
- [ ] 运维团队已培训

### 上线检查清单

- [ ] 选择低流量时段进行部署
- [ ] 通知相关人员即将上线
- [ ] 准备回滚方案
- [ ] 监控系统就绪
- [ ] 技术团队在线待命

### 上线后验证 ⚠️

- [ ] 所有服务运行正常
- [ ] 用户可以正常访问
- [ ] 关键功能测试通过
- [ ] 无严重错误日志
- [ ] 性能指标正常
- [ ] 监控告警正常

---

## 🎉 部署完成

### 完成标记

**部署日期**: _____________
**部署人员**: _____________
**生产环境URL**: _____________
**备注**: _____________

### 后续任务

- [ ] 监控系统运行24小时
- [ ] 收集用户反馈
- [ ] 性能优化调整
- [ ] 文档更新
- [ ] 团队复盘会议

---

## 📞 支持和帮助

### 文档资源

- **生产部署指南**: `PRODUCTION_DEPLOYMENT.md`
- **Phase 27报告**: `PHASE_27_COMPLETE_REPORT.md`
- **项目README**: `README.md`

### 常用命令

```bash
# 查看服务状态
./deploy-production.sh status

# 查看日志
./deploy-production.sh logs

# 重启服务
./deploy-production.sh restart

# 备份数据
./deploy-production.sh backup

# 更新部署
./deploy-production.sh update
```

### 故障处理

如遇问题，请按以下顺序排查：
1. 检查服务日志
2. 查阅故障排查文档
3. 检查配置文件
4. 联系技术负责人

---

## ⚠️ 重要提醒

1. **备份第一**: 任何重要操作前先备份
2. **测试为先**: 生产环境操作前在测试环境验证
3. **安全意识**: 定期更新密码，检查安全配置
4. **监控常态**: 保持监控系统运行，及时发现问题
5. **文档维护**: 及时更新部署文档和运维记录

---

**最后更新**: 2024-01-01
**版本**: 1.0.0
**维护者**: NGSmodule Team

✅ **检查清单使用完毕后，请保存此文件作为部署记录**
