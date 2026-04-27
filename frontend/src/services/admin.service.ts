/**
 * Admin Service - API calls for admin management
 * 对接后端 /api/v1/admin/* 端点
 */
import apiClient from './api'

// ==================== 类型定义 ====================

export interface AdminUser {
  id: string
  username: string
  email: string
  full_name?: string
  role: 'user' | 'admin'
  organization?: string
  is_active: boolean
  storage_used: number
  storage_quota: number
  storage_percent: number
  last_login?: string
  created_at: string
  updated_at: string
}

export interface AdminUserDetail extends AdminUser {
  total_projects: number
  total_samples: number
  total_tasks: number
  total_files: number
  last_activity?: string
  login_count?: number
}

export interface UserListResponse {
  users: AdminUser[]
  total: number
  page: number
  page_size: number
  total_pages: number
}

export interface UserListParams {
  skip?: number
  limit?: number
  role?: 'user' | 'admin'
  is_active?: boolean
  search?: string
  sort_by?: 'created_at' | 'username' | 'email' | 'role' | 'storage_used'
  sort_order?: 'asc' | 'desc'
}

export interface UserUpdateRequest {
  username?: string
  email?: string
  full_name?: string
  organization?: string
  storage_quota?: number
}

export interface UserActivationRequest {
  is_active: boolean
  reason?: string
}

export interface PasswordResetRequest {
  new_password: string
  notify_user?: boolean
}

export interface UserDeletionRequest {
  confirm: boolean
  transfer_data_to?: string
  reason?: string
}

export interface SystemConfig {
  general: Record<string, any>
  security: Record<string, any>
  storage: Record<string, any>
  email: Record<string, any>
  notification: Record<string, any>
  pipeline: Record<string, any>
  performance: Record<string, any>
  last_updated: string
}

export interface ConfigUpdateRequest {
  category: 'general' | 'security' | 'storage' | 'email' | 'notification' | 'pipeline' | 'performance'
  updates: Record<string, any>
  reason?: string
}

export interface LogEntry {
  timestamp: string
  level: 'debug' | 'info' | 'warning' | 'error' | 'critical'
  source: string
  message: string
  details?: Record<string, any>
  user_id?: string
  ip_address?: string
}

export interface LogResponse {
  logs: LogEntry[]
  total: number
  has_more: boolean
}

export interface ServiceHealth {
  name: string
  status: 'healthy' | 'degraded' | 'down' | 'unknown'
  response_time?: number
  last_check: string
  message?: string
  details?: Record<string, any>
}

export interface SystemHealth {
  status: 'healthy' | 'degraded' | 'down' | 'unknown'
  services: ServiceHealth[]
  timestamp: string
  uptime: number
  version: string
}

export interface CleanupOptions {
  old_logs?: boolean
  temp_files?: boolean
  failed_tasks?: boolean
  orphaned_files?: boolean
  old_notifications?: boolean
  days_to_keep?: number
  dry_run?: boolean
}

export interface CleanupResponse {
  total_items_deleted: number
  total_space_freed: number
  total_duration: number
  results: Array<{
    operation: string
    items_deleted: number
    space_freed: number
    duration: number
    errors: string[]
  }>
  timestamp: string
}

export interface AdminSystemStats {
  total_users: number
  active_users: number
  admin_users: number
  total_projects: number
  total_samples: number
  total_tasks: number
  total_files: number
  total_storage_used: number
  total_storage_allocated: number
  avg_storage_per_user: number
  tasks_today: number
  tasks_this_week: number
  tasks_this_month: number
  success_rate_today: number
  success_rate_this_week: number
  success_rate_this_month: number
  active_sessions: number
  cpu_usage: number
  memory_usage: number
  disk_usage: number
  generated_at: string
}

export interface AdminOperationResponse {
  success: boolean
  message: string
  details?: Record<string, any>
  timestamp: string
}

// ==================== Admin Service ====================

export const adminService = {
  // ============================================================
  // 用户管理
  // ============================================================

  /**
   * 获取用户列表（分页、过滤、搜索）
   */
  async getUsers(params?: UserListParams): Promise<UserListResponse> {
    return apiClient.get<UserListResponse>('/admin/users', { params })
  },

  /**
   * 获取用户详细信息
   */
  async getUser(userId: string): Promise<AdminUserDetail> {
    return apiClient.get<AdminUserDetail>(`/admin/users/${userId}`)
  },

  /**
   * 更新用户信息
   */
  async updateUser(userId: string, data: UserUpdateRequest): Promise<AdminUserDetail> {
    return apiClient.put<AdminUserDetail>(`/admin/users/${userId}`, data)
  },

  /**
   * 更改用户角色
   */
  async changeUserRole(
    userId: string,
    role: 'user' | 'admin'
  ): Promise<AdminUserDetail> {
    return apiClient.put<AdminUserDetail>(`/admin/users/${userId}/role`, { role })
  },

  /**
   * 激活/停用用户
   */
  async toggleUserStatus(
    userId: string,
    isActive: boolean,
    reason?: string
  ): Promise<AdminUserDetail> {
    return apiClient.put<AdminUserDetail>(`/admin/users/${userId}/activate`, {
      is_active: isActive,
      reason,
    } as UserActivationRequest)
  },

  /**
   * 重置用户密码
   */
  async resetUserPassword(
    userId: string,
    newPassword: string,
    notifyUser = true
  ): Promise<AdminOperationResponse> {
    return apiClient.post<AdminOperationResponse>(
      `/admin/users/${userId}/reset-password`,
      {
        new_password: newPassword,
        notify_user: notifyUser,
      } as PasswordResetRequest
    )
  },

  /**
   * 删除用户（含数据转移选项）
   */
  async deleteUser(
    userId: string,
    options: { transferDataTo?: string; reason?: string } = {}
  ): Promise<AdminOperationResponse> {
    return apiClient.delete<AdminOperationResponse>(`/admin/users/${userId}`, {
      data: {
        confirm: true,
        transfer_data_to: options.transferDataTo,
        reason: options.reason,
      } as UserDeletionRequest,
    } as any)
  },

  // ============================================================
  // 系统配置
  // ============================================================

  /**
   * 获取系统配置
   */
  async getConfig(): Promise<SystemConfig> {
    return apiClient.get<SystemConfig>('/admin/config')
  },

  /**
   * 更新系统配置
   */
  async updateConfig(request: ConfigUpdateRequest): Promise<SystemConfig> {
    return apiClient.put<SystemConfig>('/admin/config', request)
  },

  /**
   * 重置系统配置
   */
  async resetConfig(categories?: string[]): Promise<SystemConfig> {
    return apiClient.post<SystemConfig>('/admin/config/reset', {
      categories,
      confirm: true,
    })
  },

  // ============================================================
  // 系统日志
  // ============================================================

  /**
   * 查询系统日志
   */
  async getLogs(params?: {
    start_date?: string
    end_date?: string
    levels?: string[]
    sources?: string[]
    search?: string
    limit?: number
    offset?: number
  }): Promise<LogResponse> {
    return apiClient.get<LogResponse>('/admin/logs', { params })
  },

  /**
   * 下载系统日志
   */
  async downloadLogs(params?: {
    start_date?: string
    end_date?: string
    levels?: string[]
    sources?: string[]
    format?: 'json' | 'csv' | 'txt'
  }): Promise<Blob> {
    return apiClient.get('/admin/logs/download', {
      params,
      responseType: 'blob',
    } as any)
  },

  // ============================================================
  // 系统健康与维护
  // ============================================================

  /**
   * 获取系统健康状态
   */
  async getSystemHealth(): Promise<SystemHealth> {
    return apiClient.get<SystemHealth>('/admin/system/health')
  },

  /**
   * 系统清理
   */
  async cleanupSystem(options: CleanupOptions): Promise<CleanupResponse> {
    return apiClient.post<CleanupResponse>('/admin/system/cleanup', options)
  },

  /**
   * 获取系统统计（管理员）
   */
  async getSystemStats(): Promise<AdminSystemStats> {
    return apiClient.get<AdminSystemStats>('/admin/system/stats')
  },

  // ============================================================
  // 兼容旧接口（向后兼容）
  // ============================================================

  /**
   * @deprecated 使用 getUser 替代
   */
  async getUserStats(userId: string): Promise<AdminUserDetail> {
    return this.getUser(userId)
  },
}

export default adminService
