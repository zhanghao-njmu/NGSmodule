/**
 * Statistics Service
 * 统计数据API服务
 */
import apiClient from './api'

// 类型定义
export interface ProjectStats {
  total_projects: number
  active_projects: number
  archived_projects: number
  total_samples: number
  total_tasks: number
  running_tasks: number
  completed_tasks: number
  failed_tasks: number
}

export interface TaskStats {
  total_tasks: number
  pending_tasks: number
  running_tasks: number
  completed_tasks: number
  failed_tasks: number
  cancelled_tasks: number
}

export interface UserStats {
  user_id: string
  username: string
  total_projects: number
  total_samples: number
  total_tasks: number
  completed_tasks: number
  failed_tasks: number
  storage_used: number
  storage_quota: number
}

export interface SystemStats {
  total_users: number
  active_users: number
  total_projects: number
  total_samples: number
  total_tasks: number
  running_tasks: number
  completed_tasks: number
  failed_tasks: number
  total_storage_used: number
}

export interface DashboardStats {
  totalProjects: number
  activeProjects: number
  totalSamples: number
  runningTasks: number
  completedTasks: number
  failedTasks: number
  storageUsed: number
  storageQuota: number
}

class StatsService {
  /**
   * 获取项目统计
   */
  async getProjectStats(): Promise<ProjectStats> {
    const response = await apiClient.get<ProjectStats>('/items/stats')
    return response
  }

  /**
   * 获取任务统计
   */
  async getTaskStats(projectId?: string): Promise<TaskStats> {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get<TaskStats>('/tasks/stats', { params })
    return response
  }

  /**
   * 获取用户统计（需要管理员权限）
   */
  async getUserStats(userId: string): Promise<UserStats> {
    const response = await apiClient.get<UserStats>(`/users/${userId}/stats`)
    return response
  }

  /**
   * 获取系统统计（需要管理员权限）
   */
  async getSystemStats(): Promise<SystemStats> {
    const response = await apiClient.get<SystemStats>('/users/stats/system')
    return response
  }

  /**
   * 获取Dashboard所需的所有统计数据
   */
  async getDashboardStats(): Promise<DashboardStats> {
    try {
      // 并行请求项目和任务统计
      const [projectStats, _taskStats] = await Promise.all([this.getProjectStats(), this.getTaskStats()])

      return {
        totalProjects: projectStats.total_projects,
        activeProjects: projectStats.active_projects,
        totalSamples: projectStats.total_samples,
        runningTasks: projectStats.running_tasks,
        completedTasks: projectStats.completed_tasks,
        failedTasks: projectStats.failed_tasks,
        // 存储信息从用户对象中获取
        storageUsed: 0, // 将由组件从authStore获取
        storageQuota: 0, // 将由组件从authStore获取
      }
    } catch (error) {
      console.error('Failed to fetch dashboard stats:', error)
      throw error
    }
  }
}

export const statsService = new StatsService()
export default statsService
