/**
 * File Service
 * API calls for file management
 * Refactored to use CRUD factory pattern with specialized file operations
 */

import { createCrudService, extendService } from './crud.factory'
import apiClient from './api'
import type { FileItem, FileDownloadResponse, FileUploadRequest } from '@/types/file'
import type { PaginatedResponse } from '@/types/common'

/**
 * Base CRUD service for files
 * Provides getAll, getById, delete operations
 * Note: Files don't use create/update as they're managed via upload
 */
const baseCrudService = createCrudService<FileItem>({
  endpoint: 'files',
})

/**
 * Extended file service with specialized file operations
 */
export const fileService = extendService(baseCrudService, {
  /**
   * Upload a file
   * Creates a new file entry and uploads content to storage
   *
   * @param request - Upload request with sample_id, file, and progress callback
   * @returns Uploaded file item
   */
  async uploadFile({ sample_id, file, onProgress }: FileUploadRequest): Promise<FileItem> {
    const formData = new FormData()
    formData.append('file', file)

    return apiClient.post<FileItem>(`/files/upload?sample_id=${sample_id}`, formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      onUploadProgress: (progressEvent) => {
        if (onProgress && progressEvent.total) {
          const percent = Math.round((progressEvent.loaded * 100) / progressEvent.total)
          onProgress(percent)
        }
      },
    })
  },

  /**
   * Batch upload multiple files
   * Uploads multiple files for a sample
   *
   * @param sampleId - Sample ID
   * @param files - Array of files to upload
   * @param onProgress - Progress callback (optional)
   * @returns Array of uploaded file items
   */
  async batchUpload(sampleId: string, files: File[], onProgress?: (percent: number) => void): Promise<FileItem[]> {
    const formData = new FormData()
    files.forEach((file) => formData.append('files', file))

    return apiClient.post<FileItem[]>(`/files/batch-upload?sample_id=${sampleId}`, formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
      onUploadProgress: (progressEvent) => {
        if (onProgress && progressEvent.total) {
          const percent = Math.round((progressEvent.loaded * 100) / progressEvent.total)
          onProgress(percent)
        }
      },
    })
  },

  /**
   * Get file download URL
   * Returns a temporary signed URL for downloading
   *
   * @param id - File ID
   * @returns Download response with URL and filename
   */
  async getDownloadUrl(id: string): Promise<FileDownloadResponse> {
    return apiClient.get<FileDownloadResponse>(`/files/${id}/download`)
  },

  /**
   * Download file (triggers browser download)
   * Gets download URL and triggers browser download dialog
   *
   * @param id - File ID
   * @param filename - Filename for download
   */
  async downloadFile(id: string, filename: string): Promise<void> {
    const { download_url } = await this.getDownloadUrl(id)

    // Create temporary link and trigger download
    const link = document.createElement('a')
    link.href = download_url
    link.download = filename
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
  },

  /**
   * Get files by sample ID
   * Convenience method to filter files by sample
   *
   * @param sampleId - Sample ID
   * @param params - Additional query parameters
   * @returns Paginated list of files
   */
  async getFilesBySample(
    sampleId: string,
    params?: { file_type?: string; skip?: number; limit?: number },
  ): Promise<PaginatedResponse<FileItem>> {
    return baseCrudService.getAll({
      ...params,
      sample_id: sampleId,
    })
  },

  /**
   * Get files by project ID
   * Convenience method to filter files by project
   *
   * @param projectId - Project ID
   * @param params - Additional query parameters
   * @returns Paginated list of files
   */
  async getFilesByProject(
    projectId: string,
    params?: { file_type?: string; skip?: number; limit?: number },
  ): Promise<PaginatedResponse<FileItem>> {
    return baseCrudService.getAll({
      ...params,
      project_id: projectId,
    })
  },

  /**
   * Verify file checksum
   * Verifies file integrity by comparing checksums
   *
   * @param id - File ID
   * @returns Verification result
   */
  async verifyChecksum(id: string): Promise<{ valid: boolean; message: string }> {
    return apiClient.post<{ valid: boolean; message: string }>(`/files/${id}/verify`)
  },
})

export default fileService
