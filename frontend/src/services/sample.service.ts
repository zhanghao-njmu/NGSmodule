/**
 * Sample Service - API calls for sample management
 */
import apiClient from './api'
import type {
  Sample,
  SampleCreate,
  SampleUpdate,
  SampleBatchCreate,
  SampleListResponse,
} from '../types/sample'

export const sampleService = {
  /**
   * Get all samples
   */
  async getSamples(params?: {
    project_id?: string
    skip?: number
    limit?: number
  }): Promise<SampleListResponse> {
    return apiClient.get('/samples', { params })
  },

  /**
   * Get sample by ID
   */
  async getSample(id: string): Promise<Sample> {
    return apiClient.get(`/samples/${id}`)
  },

  /**
   * Create new sample
   */
  async createSample(data: SampleCreate): Promise<Sample> {
    return apiClient.post('/samples', data)
  },

  /**
   * Batch create samples
   */
  async batchCreateSamples(data: SampleBatchCreate): Promise<{ message: string; count: number }> {
    return apiClient.post('/samples/batch', data)
  },

  /**
   * Import samples from CSV
   */
  async importFromCSV(projectId: string, file: File): Promise<{ message: string; count: number }> {
    const formData = new FormData()
    formData.append('file', file)

    return apiClient.post(`/samples/import-csv?project_id=${projectId}`, formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    })
  },

  /**
   * Update sample
   */
  async updateSample(id: string, data: SampleUpdate): Promise<Sample> {
    return apiClient.put(`/samples/${id}`, data)
  },

  /**
   * Delete sample
   */
  async deleteSample(id: string): Promise<void> {
    return apiClient.delete(`/samples/${id}`)
  },
}

export default sampleService
