# Backend Integration Tests

## 概述

本目录包含NGSmodule后端API的集成测试套件，专注于验证Phase 8中实现的优化：

- **N+1查询优化**: 验证项目列表、样本列表和用户统计的查询优化
- **权限验证依赖**: 测试新的资源所有权验证依赖函数
- **全局异常处理**: 验证统一的错误响应格式
- **API响应格式**: 确保分页和响应格式的一致性

## 测试文件

### conftest.py
测试配置和fixtures:
- `db_session`: 每个测试的独立数据库会话
- `client`: FastAPI测试客户端
- `test_user`: 测试用户
- `admin_user`: 管理员用户
- `auth_headers`: 认证头
- `admin_headers`: 管理员认证头

### test_projects_api.py
项目API测试:
- ✅ 创建、读取、更新、删除项目
- ✅ N+1查询优化验证（最重要）
- ✅ 权限验证
- ✅ 错误处理
- ✅ 统计信息

### test_users_api.py
用户API测试:
- ✅ 用户信息管理
- ✅ 统计查询优化（CASE语句）
- ✅ 管理员权限
- ✅ 认证和授权

### test_exception_handlers.py
异常处理测试:
- ✅ IntegrityError处理
- ✅ ValidationError处理
- ✅ 404/401/403错误处理
- ✅ 错误响应格式一致性
- ✅ 权限依赖函数

## 安装

```bash
# 安装测试依赖
pip install -r requirements-test.txt

# 或者使用主requirements.txt（包含测试依赖）
pip install -r ../requirements.txt
```

## 运行测试

### 运行所有测试
```bash
cd backend
pytest tests/ -v
```

### 运行特定测试文件
```bash
pytest tests/test_projects_api.py -v
pytest tests/test_users_api.py -v
pytest tests/test_exception_handlers.py -v
```

### 运行特定测试类或函数
```bash
# 运行特定测试类
pytest tests/test_projects_api.py::TestProjectsAPI -v

# 运行特定测试函数
pytest tests/test_projects_api.py::TestProjectsAPI::test_list_projects_with_counts -v
```

### 生成覆盖率报告
```bash
pytest tests/ --cov=app --cov-report=html --cov-report=term
```

### 查看详细输出
```bash
pytest tests/ -v -s
```

## 重点测试场景

### 1. N+1查询优化验证

**测试**: `test_projects_api.py::TestProjectsAPI::test_list_projects_with_counts`

这个测试验证了项目列表API的N+1查询优化:
- 创建5个项目，每个有3个样本和2个任务
- 调用GET /api/v1/projects
- 验证响应包含正确的sample_count和task_count
- **关键**: 这个测试间接验证了查询数从1+2N优化到3个

### 2. 统计查询优化验证

**测试**: `test_users_api.py::TestUsersAPI::test_get_user_stats_optimized`

验证用户统计API的查询优化（从5个查询优化到3个）:
- 创建测试数据（项目、样本、不同状态的任务）
- 调用GET /api/v1/users/{id}/stats
- 验证所有统计数据正确
- **关键**: 使用CASE语句合并任务状态统计

### 3. 权限依赖函数测试

**测试**: `test_exception_handlers.py::TestPermissionDependencies`

验证新的权限验证依赖函数:
- 测试用户只能访问自己的资源
- 测试管理员可以访问所有资源
- 验证404错误（访问被拒绝）

### 4. 异常处理测试

**测试**: `test_exception_handlers.py::TestGlobalExceptionHandlers`

验证全局异常处理器:
- IntegrityError → 409 Conflict
- ValidationError → 422 Unprocessable Entity
- 一般错误 → 500 Internal Server Error
- 错误响应格式一致性

## 测试数据库

测试使用内存SQLite数据库（`:memory:`），每个测试都有独立的数据库会话:
- 测试之间完全隔离
- 无需清理
- 快速执行

## 预期结果

所有测试应该通过：
```
================================ test session starts =================================
collected 30+ items

tests/test_projects_api.py::TestProjectsAPI::test_create_project PASSED        [  3%]
tests/test_projects_api.py::TestProjectsAPI::test_list_projects_empty PASSED   [  6%]
tests/test_projects_api.py::TestProjectsAPI::test_list_projects_with_counts PASSED [10%]
...
================================ 30+ passed in 2.5s =================================
```

## 故障排除

### 导入错误
确保从backend目录运行pytest:
```bash
cd backend
pytest tests/
```

### 数据库错误
检查模型定义和关系是否正确。

### 认证错误
检查JWT配置和token生成逻辑。

## 持续集成

这些测试可以集成到CI/CD流程中：

```yaml
# .github/workflows/test.yml
- name: Run backend tests
  run: |
    cd backend
    pytest tests/ --cov=app --cov-report=xml
```

## 下一步

- [ ] 添加性能基准测试
- [ ] 添加并发测试
- [ ] 添加端到端测试
- [ ] 增加测试覆盖率到80%+

## 贡献

添加新测试时:
1. 遵循现有测试结构
2. 使用描述性的测试函数名
3. 添加文档字符串说明测试目的
4. 确保测试独立且可重复
