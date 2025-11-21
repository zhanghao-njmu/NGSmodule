/**
 * Sample types and interfaces
 */

export interface Sample {
  id: string
  project_id: string
  sample_name: string
  sample_type?: string
  metadata: Record<string, any>
  created_at: string
  updated_at: string
  file_count?: number
}

export interface SampleCreate {
  project_id: string
  sample_name: string
  sample_type?: string
  metadata?: Record<string, any>
}

export interface SampleUpdate {
  sample_name?: string
  sample_type?: string
  metadata?: Record<string, any>
}

export interface SampleBatchCreate {
  project_id: string
  samples: Array<{
    sample_name: string
    sample_type?: string
    metadata?: Record<string, any>
  }>
}

export interface SampleListResponse {
  total: number
  items: Sample[]
}
