/**
 * Project Service
 * API calls for project management
 * Refactored to use CRUD factory pattern
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { Project, ProjectCreate, ProjectUpdate, ProjectStats } from '@/types/project'

/**
 * Base CRUD service for items
 */
const baseCrudService = createCrudService<Project, ProjectCreate, ProjectUpdate>({
  endpoint: 'items',
})

/**
 * Extended project service with domain-specific methods
 */
export const projectService = extendService(baseCrudService, {
  /**
   * Get project statistics
   * Returns aggregated stats for all items
   */
  async getStats(): Promise<ProjectStats> {
    return apiClient.get<ProjectStats>('/items/stats')
  },

  /**
   * Archive a project
   * Moves project to archived status
   *
   * @param id - Project ID
   * @returns Success message
   */
  async archiveProject(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/items/${id}/archive`)
  },

  /**
   * Restore an archived project
   * Moves project back to active status
   *
   * @param id - Project ID
   * @returns Success message
   */
  async restoreProject(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/items/${id}/restore`)
  },

  /**
   * Get project with detailed information
   * Includes related samples, tasks, and statistics
   *
   * @param id - Project ID
   * @returns Project with detailed information
   */
  async getProjectDetails(id: string): Promise<
    Project & {
      samples?: any[]
      tasks?: any[]
      stats?: {
        sample_count: number
        task_count: number
        storage_used: number
      }
    }
  > {
    return apiClient.get(`/items/${id}/details`)
  },
})

export default projectService
