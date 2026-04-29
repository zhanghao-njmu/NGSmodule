/**
 * Project Service
 * API calls for project management
 * Refactored to use CRUD factory pattern
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { Project, ProjectCreate, ProjectUpdate, ProjectStats } from '@/types/project'

// The backend registers /api/v1/projects/* (see backend/app/api/v1/projects.py).
// An earlier rename of this UI page from "Projects" -> "Items" left the
// service-layer endpoint pointed at /items, which 404s in production.
// The route table for *frontend* navigation still uses /items (see App.tsx)
// but that's fine — those are React-Router paths, unrelated to API URLs.

const baseCrudService = createCrudService<Project, ProjectCreate, ProjectUpdate>({
  endpoint: 'projects',
})

export const projectService = extendService(baseCrudService, {
  async getStats(): Promise<ProjectStats> {
    return apiClient.get<ProjectStats>('/projects/stats')
  },

  async archiveProject(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/projects/${id}/archive`)
  },

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
    return apiClient.get(`/projects/${id}/details`)
  },
})

export default projectService
