/**
 * Statistics Service
 * 统计数据API服务 - 对接后端 /api/v1/stats/* 端点
 */
import apiClient from './api'

// 类型定义
export interface ProjectStats {
  total: number
  active: number
  completed: number
  failed: number
}

export interface SampleStats {
  total: number
  processed: number
  processing: number
  failed: number
}

export interface TaskStats {
  total: number
  pending: number
  running: number
  completed: number
  failed: number
  success_rate: number
}

export interface FileStats {
  total: number
  size: number
  by_type: Record<string, number>
}

export interface StorageStats {
  total: number
  used: number
  available: number
  percent_used: number
}

export interface UserStats {
  total_users: number
  active_users: number
  new_users_this_week: number
}

export interface PipelineStats {
  total_pipelines: number
  most_used: Array<{ name: string; count: number }>
}

export interface SystemStats {
  cpu_usage: number
  memory_usage: number
  disk_usage: number
  active_connections: number
  uptime_seconds: number
}

export interface StatsSummary {
  projects: ProjectStats
  samples: SampleStats
  tasks: TaskStats
  files: FileStats
  storage: StorageStats
  users?: UserStats
  pipelines?: PipelineStats
  system?: SystemStats
  generated_at: string
}

export interface QuickStats {
  total_projects: number
  active_projects: number
  total_samples: number
  total_tasks: number
  running_tasks: number
  completed_tasks: number
  storage_used: number
  storage_quota: number
  storage_percent: number
  unread_notifications: number
}

export interface TrendDataPoint {
  timestamp: string
  value: number
  label?: string
}

export interface TrendData {
  metric: string
  period: string
  data: TrendDataPoint[]
  total: number
}

// 兼容旧版的 DashboardStats 接口
export interface DashboardStats {
  totalProjects: number
  activeProjects: number
  totalSamples: number
  runningTasks: number
  completedTasks: number
  failedTasks: number
  storageUsed: number
  storageQuota: number
  unreadNotifications?: number
}

class StatsService {
  /**
   * 获取综合统计摘要
   */
  async getSummary(includeSystem = false): Promise<StatsSummary> {
    return apiClient.get<StatsSummary>('/stats/summary', {
      params: { include_system: includeSystem },
    })
  }

  /**
   * 获取项目统计
   */
  async getProjectStats(): Promise<ProjectStats> {
    return apiClient.get<ProjectStats>('/stats/projects')
  }

  /**
   * 获取样本统计
   */
  async getSampleStats(): Promise<SampleStats> {
    return apiClient.get<SampleStats>('/stats/samples')
  }

  /**
   * 获取任务统计
   */
  async getTaskStats(): Promise<TaskStats> {
    return apiClient.get<TaskStats>('/stats/tasks')
  }

  /**
   * 获取文件统计
   */
  async getFileStats(): Promise<FileStats> {
    return apiClient.get<FileStats>('/stats/files')
  }

  /**
   * 获取存储统计
   */
  async getStorageStats(): Promise<StorageStats> {
    return apiClient.get<StorageStats>('/stats/storage')
  }

  /**
   * 获取用户活动统计（需要管理员权限）
   */
  async getUserStats(): Promise<UserStats> {
    return apiClient.get<UserStats>('/stats/users')
  }

  /**
   * 获取Pipeline统计
   */
  async getPipelineStats(): Promise<PipelineStats> {
    return apiClient.get<PipelineStats>('/stats/pipelines')
  }

  /**
   * 获取系统统计（需要管理员权限）
   */
  async getSystemStats(): Promise<SystemStats> {
    return apiClient.get<SystemStats>('/stats/system')
  }

  /**
   * 获取Dashboard快速统计
   * 优化的端点，专为仪表板设计
   */
  async getQuickStats(): Promise<QuickStats> {
    return apiClient.get<QuickStats>('/stats/quick')
  }

  /**
   * 获取趋势数据
   * @param metric - 指标类型: tasks, samples, storage, projects
   * @param period - 周期: hourly, daily, weekly, monthly
   * @param days - 天数（默认30）
   */
  async getTrends(
    metric: 'tasks' | 'samples' | 'storage' | 'projects',
    period: 'hourly' | 'daily' | 'weekly' | 'monthly' = 'daily',
    days = 30
  ): Promise<TrendData> {
    return apiClient.get<TrendData>(`/stats/trends/${metric}`, {
      params: { period, days },
    })
  }

  /**
   * 获取Dashboard所需的所有统计数据
   * 兼容旧版接口
   */
  async getDashboardStats(): Promise<DashboardStats> {
    try {
      const quick = await this.getQuickStats()

      return {
        totalProjects: quick.total_projects,
        activeProjects: quick.active_projects,
        totalSamples: quick.total_samples,
        runningTasks: quick.running_tasks,
        completedTasks: quick.completed_tasks,
        failedTasks: 0, // QuickStats doesn't have failed_tasks
        storageUsed: quick.storage_used,
        storageQuota: quick.storage_quota,
        unreadNotifications: quick.unread_notifications,
      }
    } catch (error) {
      console.error('Failed to fetch dashboard stats:', error)
      throw error
    }
  }
}

export const statsService = new StatsService()
export default statsService
