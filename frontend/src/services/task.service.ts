/**
 * Task Service
 * API calls for task management
 * Refactored to use CRUD factory pattern
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { Task, TaskCreate, TaskUpdate, TaskExecuteRequest, TaskStats, TaskLogResponse } from '@/types/task'
import type { PaginatedResponse } from '@/types/common'

/**
 * Base CRUD service for tasks
 */
const baseCrudService = createCrudService<Task, TaskCreate, TaskUpdate>({
  endpoint: 'tasks',
})

/**
 * Extended task service with domain-specific methods
 */
export const taskService = extendService(baseCrudService, {
  /**
   * Get task statistics
   * Returns aggregated stats for tasks
   *
   * @param params - Optional filters (project_id)
   * @returns Task statistics
   */
  async getStats(params?: { project_id?: string }): Promise<TaskStats> {
    return apiClient.get<TaskStats>('/tasks/stats', { params })
  },

  /**
   * Execute a task
   * Starts task execution with given pipeline script
   *
   * @param id - Task ID
   * @param data - Execution request with pipeline_script and config
   * @returns Success message
   */
  async executeTask(id: string, data: TaskExecuteRequest): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/tasks/${id}/execute`, data)
  },

  /**
   * Cancel a running task
   * Stops task execution
   *
   * @param id - Task ID
   * @returns Success message
   */
  async cancelTask(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/tasks/${id}/cancel`)
  },

  /**
   * Get task execution logs
   * Retrieves logs from log file
   *
   * @param id - Task ID
   * @returns Task logs
   */
  async getTaskLogs(id: string): Promise<TaskLogResponse> {
    return apiClient.get<TaskLogResponse>(`/tasks/${id}/logs`)
  },

  /**
   * Retry a failed task
   * Re-executes a failed task with same configuration
   *
   * @param id - Task ID
   * @returns Success message
   */
  async retryTask(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/tasks/${id}/retry`)
  },

  /**
   * Get tasks by project ID
   * Convenience method to filter tasks by project
   *
   * @param projectId - Project ID
   * @param params - Additional query parameters
   * @returns Paginated list of tasks
   */
  async getTasksByProject(
    projectId: string,
    params?: { status?: string; task_type?: string },
  ): Promise<PaginatedResponse<Task>> {
    return baseCrudService.getAll({
      ...params,
      project_id: projectId,
    })
  },
})

export default taskService
