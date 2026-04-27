# 端到端测试指南 (E2E Test Guide)

**版本**: 1.0
**日期**: 2025-11-22
**Phase**: Phase 9

---

## 📖 概述

本文档提供NGSmodule系统的端到端（E2E）测试指南，涵盖从用户登录到完成完整工作流的所有步骤。

### 测试目标
- ✅ 验证前后端集成正确性
- ✅ 确保用户工作流完整无中断
- ✅ 验证数据在系统中正确流转
- ✅ 发现UI和API之间的集成问题

---

## 🚀 准备工作

### 1. 启动后端服务
```bash
cd backend
uvicorn app.main:app --reload --port 8000
```

验证:
- 访问 http://localhost:8000/docs
- API文档应该正常显示

### 2. 启动前端服务
```bash
cd frontend
npm install  # 首次运行
npm run dev
```

验证:
- 访问 http://localhost:3000 或 http://localhost:5173
- 应该看到登录页面

### 3. 准备测试数据
- 测试用户账号
- 测试CSV文件（用于样本导入）
- 测试文件（用于文件上传）

---

## 🧪 测试场景

## 场景 1: 用户认证流程

### 1.1 用户注册 (如果实现)
**前置条件**: 无

**步骤**:
1. 打开注册页面
2. 填写用户名: `testuser001`
3. 填写邮箱: `testuser001@example.com`
4. 填写密码: `Test123456!`
5. 点击"Register"按钮

**预期结果**:
- ✅ 显示成功消息
- ✅ 自动跳转到登录页或Dashboard
- ✅ 可以使用新账号登录

**检查点**:
- [ ] 用户名验证（长度、字符）
- [ ] 邮箱格式验证
- [ ] 密码强度验证
- [ ] 重复用户名提示

---

### 1.2 用户登录
**前置条件**: 用户账号已存在

**步骤**:
1. 打开登录页面 http://localhost:3000/login
2. 输入用户名: `testuser`
3. 输入密码: `testpass123`
4. 点击"Login"按钮

**预期结果**:
- ✅ 登录成功
- ✅ 跳转到Dashboard页面
- ✅ 顶部显示用户名
- ✅ Token存储在localStorage

**网络检查**:
```
POST /api/v1/auth/login
Status: 200 OK
Response:
{
  "access_token": "eyJ...",
  "token_type": "bearer"
}
```

**检查点**:
- [ ] 错误密码显示错误提示
- [ ] 不存在的用户显示错误提示
- [ ] 记住密码功能（如果有）

---

### 1.3 登出
**前置条件**: 用户已登录

**步骤**:
1. 点击顶部导航栏用户菜单
2. 点击"Logout"

**预期结果**:
- ✅ 跳转回登录页面
- ✅ localStorage清空
- ✅ 再次访问需要认证的页面会重定向到登录

---

## 场景 2: 项目管理完整流程

### 2.1 创建项目
**前置条件**: 用户已登录

**步骤**:
1. 导航到Projects页面
2. 点击"New Project"按钮
3. 填写表单:
   - Name: `E2E Test Project`
   - Description: `This is an end-to-end test project`
   - Project Type: `WGS`
4. 点击"Create"按钮

**预期结果**:
- ✅ 显示成功toast："保存成功" 或 "项目创建成功"
- ✅ 模态框关闭
- ✅ 新项目出现在列表中
- ✅ 统计卡片更新（Total +1, Active +1）

**网络检查**:
```
POST /api/v1/projects
Status: 201 Created
Response:
{
  "id": "uuid...",
  "name": "E2E Test Project",
  "description": "This is an end-to-end test project",
  "project_type": "WGS",
  "status": "active",
  "sample_count": 0,
  "task_count": 0,
  ...
}
```

**检查点**:
- [ ] 必填字段验证
- [ ] 字符长度限制
- [ ] 创建后列表正确排序（最新的在前）
- [ ] sample_count和task_count初始为0

---

### 2.2 查看项目列表
**前置条件**: 至少有1个项目

**步骤**:
1. 在Projects页面查看项目列表

**预期结果**:
- ✅ 所有项目正确显示
- ✅ 每个项目显示:
  - 项目名称（带文件夹图标）
  - 描述
  - 状态Tag
  - Sample count
  - Task count
  - 创建时间（相对时间）
  - 操作菜单（三个点）

**网络检查** (重要 - N+1优化验证):
```
GET /api/v1/projects?skip=0&limit=20
Status: 200 OK

检查Network面板:
- 应该只有3-4个SQL查询（而非1+2N）
- 响应时间应该 < 500ms（100个项目）
```

**性能检查**:
- [ ] 打开Chrome DevTools → Network
- [ ] 刷新页面
- [ ] 检查API调用次数（应该只有1个/api/v1/projects请求）
- [ ] 检查响应时间

---

### 2.3 编辑项目
**前置条件**: 项目已创建

**步骤**:
1. 点击项目行的操作菜单（三个点）
2. 选择"Edit"
3. 修改描述: `Updated description for E2E testing`
4. 点击"Save"

**预期结果**:
- ✅ 显示成功toast："更新成功"
- ✅ 模态框关闭
- ✅ 列表中的描述更新

**网络检查**:
```
PUT /api/v1/projects/{project_id}
Status: 200 OK
```

---

### 2.4 归档项目
**前置条件**: 项目已创建

**步骤**:
1. 点击项目操作菜单
2. 选择"Archive"

**预期结果**:
- ✅ 显示loading toast："归档中..."
- ✅ 成功后显示："项目已归档"
- ✅ 项目状态变为"archived"
- ✅ 项目名称旁显示"archived" tag
- ✅ 统计卡片更新（Active -1, Archived +1）

**网络检查**:
```
PUT /api/v1/projects/{project_id}/archive
Status: 200 OK
```

---

### 2.5 恢复项目
**前置条件**: 项目已归档

**步骤**:
1. 找到归档的项目
2. 点击操作菜单
3. 选择"Restore"

**预期结果**:
- ✅ 显示loading toast："恢复中..."
- ✅ 成功后显示："项目已恢复"
- ✅ 项目状态变回"active"
- ✅ "archived" tag消失
- ✅ 统计卡片更新

---

### 2.6 删除项目 ⚠️ 重要
**前置条件**: 项目已创建

**步骤**:
1. 点击项目操作菜单
2. 选择"Delete"
3. **验证ConfirmDialog显示**:
   - 标题应该是"删除项目"或类似文本
   - 内容应包含项目名称
   - 应有警告信息
4. 点击"取消"（第一次）
5. 再次打开删除对话框
6. 点击"删除"（第二次）

**预期结果**:
- ✅ 第一次点击"取消"：对话框关闭，项目仍在列表中
- ✅ 第二次点击"删除"：
  - 显示loading toast："删除中..."
  - 成功后显示："删除成功"
  - 项目从列表中消失
  - 统计卡片更新（Total -1）

**网络检查**:
```
DELETE /api/v1/projects/{project_id}
Status: 200 or 204
```

**检查点**:
- [ ] ConfirmDialog样式正确（危险操作应该是红色按钮）
- [ ] 删除是真删除还是软删除
- [ ] 如果项目有关联数据（样本、任务），应该一并删除或阻止删除

---

## 场景 3: 样本管理流程

### 3.1 导航到样本页面
**前置条件**: 用户已登录

**步骤**:
1. 点击侧边栏"Samples"菜单
2. 观察页面加载

**预期结果**:
- ✅ 页面正常加载
- ✅ 显示项目选择下拉框
- ✅ 显示"Select a project"提示
- ✅ 显示空状态（未选择项目时）

---

### 3.2 选择项目
**前置条件**: 至少有1个项目

**步骤**:
1. 点击项目选择下拉框
2. 选择"E2E Test Project"

**预期结果**:
- ✅ 显示loading状态
- ✅ 加载该项目的样本列表
- ✅ 如果没有样本，显示空状态提示

**网络检查**:
```
GET /api/v1/samples?project_id={project_id}
Status: 200 OK
```

---

### 3.3 CSV导入样本 ⚠️ 重要
**前置条件**: 已选择项目，准备好CSV文件

**测试CSV内容** (samples.csv):
```csv
sample_name,sample_type,group_name
Sample001,WGS,Control
Sample002,WGS,Treatment
Sample003,RNA-Seq,Control
```

**步骤**:
1. 点击"Import CSV"按钮
2. 选择准备好的CSV文件
3. 观察导入过程

**预期结果**:
- ✅ 点击后立即显示loading toast："导入中..."
- ✅ 导入成功后:
  - Loading toast消失
  - 显示成功toast："上传成功"
  - 样本列表刷新
  - 新样本出现在列表中（3个样本）
- ✅ 导入失败时:
  - 显示错误toast："上传失败，请重试"

**网络检查**:
```
POST /api/v1/samples/import
Content-Type: multipart/form-data
Status: 200 OK
Response:
{
  "imported": 3,
  "failed": 0,
  "message": "Successfully imported 3 samples"
}
```

**检查点**:
- [ ] CSV格式验证
- [ ] 重复样本名称处理
- [ ] 必填字段验证
- [ ] 错误行报告

---

### 3.4 查看样本列表
**前置条件**: 项目中有样本

**步骤**:
1. 查看样本列表

**预期结果**:
- ✅ 所有样本正确显示
- ✅ 每个样本显示:
  - 样本名称（带实验图标）
  - 样本类型（Tag）
  - File count
  - 创建时间
- ✅ file_count正确（初始为0）

**网络检查** (N+1优化验证):
```
GET /api/v1/samples?project_id={project_id}
应该只有2个SQL查询（不是1+N）
```

---

## 场景 4: 任务管理流程

### 4.1 查看任务列表
**前置条件**: 用户已登录

**步骤**:
1. 点击侧边栏"Tasks"菜单
2. 观察页面加载

**预期结果**:
- ✅ 页面正常加载
- ✅ 统计卡片显示（Total, Running, Completed, Failed）
- ✅ 任务列表显示
- ✅ 控制台显示WebSocket连接成功

**WebSocket检查**:
```
开发者工具 → Network → WS
应该看到WebSocket连接
Status: 101 Switching Protocols
```

---

### 4.2 创建任务
**前置条件**: 有项目和样本

**步骤**:
1. 点击"New Task"按钮（如果有）
2. 填写任务信息
3. 提交创建

**预期结果**:
- ✅ 任务创建成功
- ✅ 出现在任务列表中
- ✅ 状态为"pending"或"running"

---

### 4.3 监控任务执行
**前置条件**: 有运行中的任务

**步骤**:
1. 观察运行中的任务
2. 保持页面打开（不刷新）

**预期结果**:
- ✅ 任务状态自动更新（通过WebSocket）
- ✅ 进度条实时更新
- ✅ 统计卡片实时更新
- ✅ 无需手动刷新页面

**WebSocket消息检查**:
```
应该收到任务更新消息:
{
  "type": "task_update",
  "task_id": "...",
  "status": "running",
  "progress": 45
}
```

---

### 4.4 取消任务 ⚠️ 重要
**前置条件**: 有运行中的任务

**步骤**:
1. 找到状态为"running"的任务
2. 点击"Cancel"按钮
3. **验证Confirm对话框**:
   - 标题："取消任务"
   - 内容包含任务名称
4. 点击"保持运行"（第一次）
5. 再次点击Cancel
6. 点击"取消任务"（第二次）

**预期结果**:
- ✅ 第一次点"保持运行"：对话框关闭，任务继续运行
- ✅ 第二次点"取消任务"：
  - 显示loading toast："取消中..."
  - 成功后显示："任务已取消"
  - 任务状态变为"cancelled"
  - 统计卡片更新（Running -1）

**网络检查**:
```
PUT /api/v1/tasks/{task_id}/cancel
Status: 200 OK
```

---

## 场景 5: 权限和安全测试

### 5.1 访问其他用户的资源
**前置条件**:
- 创建两个用户账号
- User A创建了一个项目

**步骤**:
1. 用User A登录，创建项目，记下项目ID
2. 登出
3. 用User B登录
4. 尝试通过直接URL访问User A的项目
   - URL: `/projects/{user_a_project_id}`

**预期结果**:
- ✅ 应该看到404错误或"Access Denied"
- ✅ User B不能看到User A的项目
- ✅ User B不能编辑User A的项目

**网络检查**:
```
GET /api/v1/projects/{user_a_project_id}
Status: 404 Not Found
Response:
{
  "detail": "Project not found or access denied"
}
```

---

### 5.2 管理员权限
**前置条件**: 有管理员账号

**步骤**:
1. 用管理员账号登录
2. 访问"Users"页面（仅管理员可见）

**预期结果**:
- ✅ 管理员可以看到所有用户
- ✅ 管理员可以看到用户统计
- ✅ 管理员可以访问所有项目

---

### 5.3 Token过期处理
**前置条件**: 用户已登录

**步骤**:
1. 登录系统
2. 在Chrome DevTools → Application → LocalStorage
3. 找到auth_token，记下它
4. 等待Token过期（或手动修改为无效token）
5. 尝试执行需要认证的操作

**预期结果**:
- ✅ 显示"您的登录已过期，请重新登录"
- ✅ 自动跳转到登录页面
- ✅ LocalStorage中的token被清除

**网络检查**:
```
任何API请求
Status: 401 Unauthorized
Response:
{
  "detail": "Could not validate credentials"
}
```

---

## 场景 6: 错误处理和恢复

### 6.1 网络错误处理
**步骤**:
1. 登录系统
2. 打开Chrome DevTools → Network
3. 设置"Offline"模式（模拟断网）
4. 尝试创建项目

**预期结果**:
- ✅ 显示错误toast："网络连接失败，请检查网络"
- ✅ 不会有无限loading
- ✅ 可以重试操作

**恢复步骤**:
1. 关闭"Offline"模式
2. 重试创建项目
3. 应该成功

---

### 6.2 服务器错误处理
**步骤**:
1. 停止后端服务器
2. 在前端尝试操作

**预期结果**:
- ✅ 显示友好的错误消息
- ✅ 不会崩溃或白屏
- ✅ 可以导航到其他页面

---

### 6.3 表单验证
**步骤**:
1. 打开创建项目表单
2. 不填写必填字段
3. 点击提交

**预期结果**:
- ✅ 显示字段验证错误
- ✅ 错误消息清晰
- ✅ 焦点自动移到第一个错误字段

---

## 📊 性能测试

### 加载时间测试
使用Chrome DevTools → Performance

1. **首页加载**:
   - 目标: < 2s
   - 记录: _____s

2. **项目列表加载** (100个项目):
   - 目标: < 500ms
   - 记录: _____ms

3. **样本列表加载** (200个样本):
   - 目标: < 500ms
   - 记录: _____ms

### SQL查询优化验证

**检查方法**:
1. 启用后端SQL日志
2. 访问项目列表页面
3. 统计SQL查询数量

**预期结果**:
- 项目列表: 3个查询（不是1+2N）
- 样本列表: 2个查询（不是1+N）
- 用户统计: 3个查询（不是5个）

---

## ✅ 测试完成清单

### 功能测试
- [ ] 用户认证流程（登录、登出）
- [ ] 项目CRUD操作
- [ ] 样本CSV导入
- [ ] 任务监控和取消
- [ ] 权限验证
- [ ] 错误处理

### UI/UX测试
- [ ] ConfirmDialog正常工作
- [ ] Toast通知正常工作
- [ ] Loading状态正确
- [ ] WebSocket实时更新

### 性能测试
- [ ] N+1查询已优化
- [ ] API响应时间达标
- [ ] 页面加载速度正常

### 安全测试
- [ ] 权限隔离正确
- [ ] Token过期处理正确
- [ ] XSS防护（输入验证）

---

## 🐛 问题报告模板

```markdown
## Bug Report

**发现时间**: [日期时间]
**场景**: [场景编号和名称]
**严重程度**: [Critical/High/Medium/Low]

**步骤重现**:
1.
2.
3.

**预期结果**:
-

**实际结果**:
-

**截图/录屏**:
[附上截图]

**浏览器信息**:
- 浏览器:
- 版本:

**控制台错误**:
```
[粘贴console错误]
```

**网络请求**:
```
[粘贴相关API请求/响应]
```
```

---

**文档版本**: 1.0
**最后更新**: 2025-11-22
