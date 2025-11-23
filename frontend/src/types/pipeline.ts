/**
 * Pipeline types and interfaces
 */

import type { PaginatedResponse } from './common'

export interface PipelineTemplate {
  id: string
  name: string
  display_name: string
  description?: string
  category: string
  script_name: string
  script_path?: string
  default_params: Record<string, any>
  param_schema: Record<string, ParamSchema>
  estimated_time?: string
  min_memory_gb?: string
  min_cpu_cores?: string
  is_active: boolean
  is_builtin: boolean
  sort_order: string
  tags: string[]
}

export interface ParamSchema {
  type: 'integer' | 'number' | 'string' | 'boolean' | 'select'
  label: string
  min?: number
  max?: number
  step?: number
  default: any
  options?: Array<string | { value: any; label: string }>
}

/**
 * @deprecated Use PaginatedResponse<PipelineTemplate> instead
 */
export type PipelineTemplateListResponse = PaginatedResponse<PipelineTemplate>

export interface PipelineTemplateCategory {
  category: string
  count: number
  templates: string[]
}

export interface PipelineExecuteRequest {
  template_id: string
  task_name: string
  project_id: string
  sample_ids: string[]
  parameters: Record<string, any>
}

export interface PipelineBatchExecuteRequest {
  template_id: string
  project_id: string
  sample_ids: string[]
  task_name_prefix: string
  parameters: Record<string, any>
}

export interface PipelineBatchExecuteResponse {
  total_tasks: number
  created_tasks: string[]
  failed_samples: Array<{
    sample_id: string
    sample_name: string
    error: string
  }>
}

export interface ParameterRecommendationResponse {
  recommended_params: Record<string, any>
  confidence_score: number
  based_on_tasks: number
  explanation: string
}
