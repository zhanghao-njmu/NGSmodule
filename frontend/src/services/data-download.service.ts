/**
 * Data downloads + vendor credentials API client.
 *
 * apiClient.{get,post,...} returns Promise<T> directly (the response
 * data, not an AxiosResponse wrapper).
 */
import apiClient from './api'
import type {
  SessionStatus,
  SessionLoginRequest,
  DownloadJob,
  DownloadJobCreateRequest,
  VendorCredential,
  VendorCredentialCreateRequest,
} from '@/types/data-download'

interface ListResponse<T> {
  total: number
  items: T[]
}

interface MessageResponse {
  message: string
}

export const dataDownloadService = {
  // Vendor catalog
  listVendors: () => apiClient.get<string[]>('/data-downloads/vendors'),

  // Sessions
  getSession: (vendor: string) => apiClient.get<SessionStatus>(`/data-downloads/sessions/${vendor}`),
  openSession: (req: SessionLoginRequest) => apiClient.post<SessionStatus>('/data-downloads/sessions', req),
  closeSession: (vendor: string) => apiClient.delete<MessageResponse>(`/data-downloads/sessions/${vendor}`),

  // Jobs
  createJob: (req: DownloadJobCreateRequest) => apiClient.post<DownloadJob>('/data-downloads/jobs', req),
  listJobs: () => apiClient.get<ListResponse<DownloadJob>>('/data-downloads/jobs'),
  getJob: (id: string) => apiClient.get<DownloadJob>(`/data-downloads/jobs/${id}`),
  cancelJob: (id: string) => apiClient.delete<DownloadJob>(`/data-downloads/jobs/${id}`),
}

export const vendorCredentialService = {
  list: (vendor?: string) =>
    apiClient.get<ListResponse<VendorCredential>>(
      vendor ? `/vendor-credentials?vendor=${encodeURIComponent(vendor)}` : '/vendor-credentials',
    ),
  create: (req: VendorCredentialCreateRequest) => apiClient.post<VendorCredential>('/vendor-credentials', req),
  delete: (id: string) => apiClient.delete<MessageResponse>(`/vendor-credentials/${id}`),
}
