/**
 * Sample Service
 * API calls for sample management
 * Refactored to use CRUD factory pattern
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { Sample, SampleCreate, SampleUpdate, SampleBatchCreate } from '@/types/sample'
import type { PaginatedResponse } from '@/types/common'

/**
 * Base CRUD service for samples
 */
const baseCrudService = createCrudService<Sample, SampleCreate, SampleUpdate>({
  endpoint: 'samples',
})

/**
 * Extended sample service with domain-specific methods
 */
export const sampleService = extendService(baseCrudService, {
  /**
   * Batch create multiple samples
   * Creates multiple samples in a single transaction
   *
   * @param data - Batch create data with project_id and array of samples
   * @returns Success message and count of created samples
   */
  async batchCreateSamples(data: SampleBatchCreate): Promise<{ message: string; count: number }> {
    return apiClient.post<{ message: string; count: number }>('/samples/batch', data)
  },

  /**
   * Import samples from CSV file
   * Parses CSV and creates samples in bulk
   *
   * @param projectId - Project ID to associate samples with
   * @param file - CSV file containing sample data
   * @returns Success message and count of imported samples
   */
  async importFromCSV(projectId: string, file: File): Promise<{ message: string; count: number }> {
    const formData = new FormData()
    formData.append('file', file)

    return apiClient.post<{ message: string; count: number }>(
      `/samples/import-csv?project_id=${projectId}`,
      formData,
      {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      }
    )
  },

  /**
   * Export samples to CSV file
   * Downloads samples for a project as CSV
   *
   * @param projectId - Project ID to export samples from
   * @returns Download URL for CSV file
   */
  async exportToCSV(projectId: string): Promise<{ download_url: string }> {
    return apiClient.get<{ download_url: string }>(`/samples/export-csv`, {
      params: { project_id: projectId },
    })
  },

  /**
   * Get samples by project ID
   * Convenience method to filter samples by project
   *
   * @param projectId - Project ID
   * @param params - Additional query parameters
   * @returns Paginated list of samples
   */
  async getSamplesByProject(
    projectId: string,
    params?: { skip?: number; limit?: number }
  ): Promise<PaginatedResponse<Sample>> {
    return baseCrudService.getAll({
      ...params,
      project_id: projectId,
    })
  },

  // Backward compatibility aliases
  // These can be removed after all components are migrated

  /**
   * @deprecated Use getAll instead
   */
  getSamples(params?: { project_id?: string; skip?: number; limit?: number }): Promise<PaginatedResponse<Sample>> {
    return baseCrudService.getAll(params)
  },

  /**
   * @deprecated Use getById instead
   */
  getSample(id: string): Promise<Sample> {
    return baseCrudService.getById(id)
  },

  /**
   * @deprecated Use create instead
   */
  createSample(data: SampleCreate): Promise<Sample> {
    return baseCrudService.create(data)
  },

  /**
   * @deprecated Use update instead
   */
  updateSample(id: string, data: SampleUpdate): Promise<Sample> {
    return baseCrudService.update(id, data)
  },

  /**
   * @deprecated Use delete instead
   */
  deleteSample(id: string): Promise<void> {
    return baseCrudService.delete(id)
  },
})

export default sampleService
