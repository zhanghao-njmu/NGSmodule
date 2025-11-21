/**
 * Task Service - API calls for task management
 */
import apiClient from './api'
import type {
  Task,
  TaskCreate,
  TaskUpdate,
  TaskExecuteRequest,
  TaskListResponse,
  TaskStats,
  TaskLogResponse,
} from '../types/task'

export const taskService = {
  /**
   * Get all tasks
   */
  async getTasks(params?: {
    project_id?: string
    status?: string
    task_type?: string
  }): Promise<TaskListResponse> {
    return apiClient.get('/tasks', { params })
  },

  /**
   * Get task statistics
   */
  async getStats(params?: { project_id?: string }): Promise<TaskStats> {
    return apiClient.get('/tasks/stats', { params })
  },

  /**
   * Get task by ID
   */
  async getTask(id: string): Promise<Task> {
    return apiClient.get(`/tasks/${id}`)
  },

  /**
   * Create new task
   */
  async createTask(data: TaskCreate): Promise<Task> {
    return apiClient.post('/tasks', data)
  },

  /**
   * Update task
   */
  async updateTask(id: string, data: TaskUpdate): Promise<Task> {
    return apiClient.put(`/tasks/${id}`, data)
  },

  /**
   * Execute task
   */
  async executeTask(id: string, data: TaskExecuteRequest): Promise<{ message: string }> {
    return apiClient.post(`/tasks/${id}/execute`, data)
  },

  /**
   * Cancel task
   */
  async cancelTask(id: string): Promise<{ message: string }> {
    return apiClient.post(`/tasks/${id}/cancel`)
  },

  /**
   * Get task logs
   */
  async getTaskLogs(id: string): Promise<TaskLogResponse> {
    return apiClient.get(`/tasks/${id}/logs`)
  },

  /**
   * Delete task
   */
  async deleteTask(id: string): Promise<void> {
    return apiClient.delete(`/tasks/${id}`)
  },
}

export default taskService
