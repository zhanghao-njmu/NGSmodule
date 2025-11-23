/**
 * Project Service
 * API calls for project management
 * Refactored to use CRUD factory pattern
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { Project, ProjectCreate, ProjectUpdate, ProjectStats } from '@/types/project'
import type { PaginatedResponse } from '@/types/common'

/**
 * Base CRUD service for projects
 */
const baseCrudService = createCrudService<Project, ProjectCreate, ProjectUpdate>({
  endpoint: 'projects',
})

/**
 * Extended project service with domain-specific methods
 */
export const projectService = extendService(baseCrudService, {
  /**
   * Get project statistics
   * Returns aggregated stats for all projects
   */
  async getStats(): Promise<ProjectStats> {
    return apiClient.get<ProjectStats>('/projects/stats')
  },

  /**
   * Archive a project
   * Moves project to archived status
   *
   * @param id - Project ID
   * @returns Success message
   */
  async archiveProject(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/projects/${id}/archive`)
  },

  /**
   * Restore an archived project
   * Moves project back to active status
   *
   * @param id - Project ID
   * @returns Success message
   */
  async restoreProject(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/projects/${id}/restore`)
  },

  /**
   * Get project with detailed information
   * Includes related samples, tasks, and statistics
   *
   * @param id - Project ID
   * @returns Project with detailed information
   */
  async getProjectDetails(id: string): Promise<Project & {
    samples?: any[]
    tasks?: any[]
    stats?: {
      sample_count: number
      task_count: number
      storage_used: number
    }
  }> {
    return apiClient.get(`/projects/${id}/details`)
  },

  // Backward compatibility aliases
  // These can be removed after all components are migrated

  /**
   * @deprecated Use getAll instead
   */
  getProjects(params?: { status?: string; skip?: number; limit?: number }): Promise<PaginatedResponse<Project>> {
    return baseCrudService.getAll(params)
  },

  /**
   * @deprecated Use getById instead
   */
  getProject(id: string): Promise<Project> {
    return baseCrudService.getById(id)
  },

  /**
   * @deprecated Use create instead
   */
  createProject(data: ProjectCreate): Promise<Project> {
    return baseCrudService.create(data)
  },

  /**
   * @deprecated Use update instead
   */
  updateProject(id: string, data: ProjectUpdate): Promise<Project> {
    return baseCrudService.update(id, data)
  },

  /**
   * @deprecated Use delete instead
   */
  deleteProject(id: string): Promise<void> {
    return baseCrudService.delete(id)
  },
})

export default projectService
