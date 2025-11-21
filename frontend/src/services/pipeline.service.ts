/**
 * Pipeline Service - API calls for pipeline management
 */
import apiClient from './api'
import type {
  PipelineTemplate,
  PipelineTemplateListResponse,
  PipelineTemplateCategory,
  PipelineExecuteRequest,
  PipelineBatchExecuteRequest,
  PipelineBatchExecuteResponse,
  ParameterRecommendationResponse,
} from '../types/pipeline'
import type { Task } from '../types/task'

export const pipelineService = {
  /**
   * Get all pipeline templates
   */
  async getTemplates(params?: {
    category?: string
    is_active?: boolean
    search?: string
  }): Promise<PipelineTemplateListResponse> {
    return apiClient.get('/pipelines', { params })
  },

  /**
   * Get pipeline categories
   */
  async getCategories(): Promise<PipelineTemplateCategory[]> {
    return apiClient.get('/pipelines/categories')
  },

  /**
   * Get pipeline template by ID
   */
  async getTemplate(id: string): Promise<PipelineTemplate> {
    return apiClient.get(`/pipelines/${id}`)
  },

  /**
   * Execute pipeline
   */
  async executePipeline(data: PipelineExecuteRequest): Promise<Task> {
    return apiClient.post('/pipelines/execute', data)
  },

  /**
   * Batch execute pipeline on multiple samples
   */
  async batchExecutePipeline(data: PipelineBatchExecuteRequest): Promise<PipelineBatchExecuteResponse> {
    return apiClient.post('/pipelines/batch-execute', data)
  },

  /**
   * Get AI-powered parameter recommendations
   */
  async getParameterRecommendations(
    templateId: string,
    projectId?: string
  ): Promise<ParameterRecommendationResponse> {
    const params = projectId ? { project_id: projectId } : {}
    return apiClient.get(`/pipelines/${templateId}/recommend-parameters`, { params })
  },

  /**
   * Toggle pipeline template active status (Admin only)
   */
  async toggleTemplate(id: string): Promise<{ message: string }> {
    return apiClient.post(`/pipelines/${id}/toggle`)
  },
}

export default pipelineService
