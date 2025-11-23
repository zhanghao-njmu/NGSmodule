/**
 * Pipeline Service
 * API calls for pipeline management and execution
 * Refactored with standardized methods
 */

import apiClient from './api'
import type {
  PipelineTemplate,
  PipelineTemplateCategory,
  PipelineExecuteRequest,
  PipelineBatchExecuteRequest,
  PipelineBatchExecuteResponse,
  ParameterRecommendationResponse,
} from '@/types/pipeline'
import type { Task } from '@/types/task'
import type { PaginatedResponse } from '@/types/common'

/**
 * Pipeline service
 * Handles pipeline templates and execution
 */
export const pipelineService = {
  /**
   * Get all pipeline templates with optional filters
   *
   * @param params - Query parameters
   * @returns Paginated list of pipeline templates
   */
  async getAll(params?: {
    category?: string
    is_active?: boolean
    search?: string
  }): Promise<PaginatedResponse<PipelineTemplate>> {
    return apiClient.get<PaginatedResponse<PipelineTemplate>>('/pipelines', { params })
  },

  /**
   * Get pipeline template by ID
   *
   * @param id - Pipeline template ID
   * @returns Pipeline template details
   */
  async getById(id: string): Promise<PipelineTemplate> {
    return apiClient.get<PipelineTemplate>(`/pipelines/${id}`)
  },

  /**
   * Get pipeline categories
   * Returns available categories with counts
   *
   * @returns Array of pipeline categories
   */
  async getCategories(): Promise<PipelineTemplateCategory[]> {
    return apiClient.get<PipelineTemplateCategory[]>('/pipelines/categories')
  },

  /**
   * Execute a pipeline
   * Creates and starts a task for pipeline execution
   *
   * @param data - Execution request with template_id, sample_id, parameters
   * @returns Created task
   */
  async executePipeline(data: PipelineExecuteRequest): Promise<Task> {
    return apiClient.post<Task>('/pipelines/execute', data)
  },

  /**
   * Batch execute pipeline on multiple samples
   * Creates multiple tasks for batch execution
   *
   * @param data - Batch execution request with template_id and sample_ids
   * @returns Batch execution response with created task IDs
   */
  async batchExecutePipeline(data: PipelineBatchExecuteRequest): Promise<PipelineBatchExecuteResponse> {
    return apiClient.post<PipelineBatchExecuteResponse>('/pipelines/batch-execute', data)
  },

  /**
   * Get AI-powered parameter recommendations
   * Returns recommended parameters based on project history and best practices
   *
   * @param templateId - Pipeline template ID
   * @param projectId - Optional project ID for context
   * @returns Parameter recommendations
   */
  async getParameterRecommendations(
    templateId: string,
    projectId?: string
  ): Promise<ParameterRecommendationResponse> {
    const params = projectId ? { project_id: projectId } : {}
    return apiClient.get<ParameterRecommendationResponse>(
      `/pipelines/${templateId}/recommend-parameters`,
      { params }
    )
  },

  /**
   * Toggle pipeline template active status
   * Admin only - enables/disables a pipeline template
   *
   * @param id - Pipeline template ID
   * @returns Success message
   */
  async toggleTemplate(id: string): Promise<{ message: string }> {
    return apiClient.post<{ message: string }>(`/pipelines/${id}/toggle`)
  },

  /**
   * Get pipelines by category
   * Convenience method to filter by category
   *
   * @param category - Pipeline category
   * @returns Paginated list of pipeline templates
   */
  async getByCategory(category: string): Promise<PaginatedResponse<PipelineTemplate>> {
    return this.getAll({ category })
  },

  /**
   * Get active pipelines only
   * Convenience method to filter active templates
   *
   * @returns Paginated list of active pipeline templates
   */
  async getActive(): Promise<PaginatedResponse<PipelineTemplate>> {
    return this.getAll({ is_active: true })
  },

  /**
   * Search pipelines
   * Convenience method for searching templates
   *
   * @param query - Search query
   * @returns Paginated list of matching pipeline templates
   */
  async search(query: string): Promise<PaginatedResponse<PipelineTemplate>> {
    return this.getAll({ search: query })
  },

  // Backward compatibility aliases

  /**
   * @deprecated Use getAll instead
   */
  getTemplates(params?: {
    category?: string
    is_active?: boolean
    search?: string
  }): Promise<PaginatedResponse<PipelineTemplate>> {
    return this.getAll(params)
  },

  /**
   * @deprecated Use getById instead
   */
  getTemplate(id: string): Promise<PipelineTemplate> {
    return this.getById(id)
  },
}

export default pipelineService
