/**
 * Type definitions for the data-downloads + vendor-credentials API.
 */

export type DownloadStatus = 'pending' | 'running' | 'completed' | 'failed' | 'cancelled'

export interface SessionStatus {
  vendor: string
  active: boolean
  pid?: number | null
}

export interface SessionLoginRequest {
  vendor: string
  email?: string
  password?: string
  credential_id?: string
}

export interface DownloadJob {
  id: string
  vendor: string
  source_path: string
  dest_path: string
  log_path?: string | null
  status: DownloadStatus
  progress_pct: number
  bytes_downloaded?: number | null
  file_size?: number | null
  error_message?: string | null
  started_at?: string | null
  finished_at?: string | null
  created_at: string
  project_id?: string | null
}

export interface DownloadJobCreateRequest {
  vendor: string
  source_path: string
  dest_path: string
  auto_register?: boolean
  project_name?: string
}

export interface VendorCredential {
  id: string
  vendor: string
  label: string
  email_preview: string
  created_at: string
  last_used_at?: string | null
}

export interface VendorCredentialCreateRequest {
  vendor: string
  label?: string
  email: string
  password: string
}
