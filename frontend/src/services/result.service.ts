/**
 * Results Service - API calls for analysis results
 */
import api from './api'
import type { Result, ResultListResponse, ResultVisualizationData, ResultSummary } from '../types/result'

class ResultService {
  private baseURL = '/results'

  /**
   * List results with filtering and pagination
   */
  async listResults(params?: {
    task_id?: string
    result_type?: string
    skip?: number
    limit?: number
  }): Promise<ResultListResponse> {
    const response = await api.get<ResultListResponse>(this.baseURL, { params })
    return response.data
  }

  /**
   * Get a specific result by ID
   */
  async getResult(resultId: string): Promise<Result> {
    const response = await api.get<Result>(`${this.baseURL}/${resultId}`)
    return response.data
  }

  /**
   * Get visualization data for a result
   */
  async getVisualizationData(resultId: string): Promise<ResultVisualizationData> {
    const response = await api.get<ResultVisualizationData>(
      `${this.baseURL}/${resultId}/visualization`
    )
    return response.data
  }

  /**
   * Get summary of all results for a task
   */
  async getTaskResultsSummary(taskId: string): Promise<ResultSummary> {
    const response = await api.get<ResultSummary>(`${this.baseURL}/task/${taskId}/summary`)
    return response.data
  }

  /**
   * Download result file
   */
  async downloadResult(result: Result): Promise<void> {
    if (!result.result_path) {
      throw new Error('Result has no downloadable file')
    }

    // In a real implementation, this would trigger file download
    window.open(result.result_path, '_blank')
  }
}

export default new ResultService()
