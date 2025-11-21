/**
 * Project Service - API calls for project management
 */
import apiClient from './api'
import type {
  Project,
  ProjectCreate,
  ProjectUpdate,
  ProjectListResponse,
  ProjectStats,
} from '../types/project'

export const projectService = {
  /**
   * Get all projects
   */
  async getProjects(params?: {
    status?: string
    skip?: number
    limit?: number
  }): Promise<ProjectListResponse> {
    return apiClient.get('/projects', { params })
  },

  /**
   * Get project statistics
   */
  async getStats(): Promise<ProjectStats> {
    return apiClient.get('/projects/stats')
  },

  /**
   * Get project by ID
   */
  async getProject(id: string): Promise<Project> {
    return apiClient.get(`/projects/${id}`)
  },

  /**
   * Create new project
   */
  async createProject(data: ProjectCreate): Promise<Project> {
    return apiClient.post('/projects', data)
  },

  /**
   * Update project
   */
  async updateProject(id: string, data: ProjectUpdate): Promise<Project> {
    return apiClient.put(`/projects/${id}`, data)
  },

  /**
   * Delete project
   */
  async deleteProject(id: string): Promise<void> {
    return apiClient.delete(`/projects/${id}`)
  },

  /**
   * Archive project
   */
  async archiveProject(id: string): Promise<{ message: string }> {
    return apiClient.post(`/projects/${id}/archive`)
  },

  /**
   * Restore archived project
   */
  async restoreProject(id: string): Promise<{ message: string }> {
    return apiClient.post(`/projects/${id}/restore`)
  },
}

export default projectService
