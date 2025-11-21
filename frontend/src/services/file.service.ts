/**
 * File Service - API calls for file management
 */
import apiClient from './api'
import type { FileItem, FileListResponse, FileDownloadResponse, FileUploadRequest } from '../types/file'

export const fileService = {
  /**
   * Get all files
   */
  async getFiles(params?: {
    sample_id?: string
    project_id?: string
    file_type?: string
    skip?: number
    limit?: number
  }): Promise<FileListResponse> {
    return apiClient.get('/files', { params })
  },

  /**
   * Get file by ID
   */
  async getFile(id: string): Promise<FileItem> {
    return apiClient.get(`/files/${id}`)
  },

  /**
   * Upload file
   */
  async uploadFile({ sample_id, file, onProgress }: FileUploadRequest): Promise<FileItem> {
    const formData = new FormData()
    formData.append('file', file)

    return apiClient.post(`/files/upload?sample_id=${sample_id}`, formData, {
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
   */
  async getDownloadUrl(id: string): Promise<FileDownloadResponse> {
    return apiClient.get(`/files/${id}/download`)
  },

  /**
   * Delete file
   */
  async deleteFile(id: string): Promise<void> {
    return apiClient.delete(`/files/${id}`)
  },

  /**
   * Download file (triggers browser download)
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
}

export default fileService
