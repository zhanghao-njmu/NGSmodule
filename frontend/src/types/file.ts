/**
 * File types and interfaces
 */

import type { PaginatedResponse } from './common'

export interface FileItem {
  id: string
  sample_id: string
  filename: string
  file_type?: string
  file_size: number
  file_path: string
  checksum?: string
  uploaded_at: string
  project_id?: string
}

export interface FileUploadRequest {
  sample_id: string
  file: File
  onProgress?: (percent: number) => void
}

/**
 * @deprecated Use PaginatedResponse<FileItem> instead
 */
export type FileListResponse = PaginatedResponse<FileItem>

export interface FileDownloadResponse {
  download_url: string
  filename: string
}
