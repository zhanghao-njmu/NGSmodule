/**
 * Sample types and interfaces
 * 与后端API保持一致
 */

import type { PaginatedResponse } from './common'

export interface Sample {
  id: string
  project_id: string
  sample_id: string // 样本标识符
  run_id?: string // 测序批次ID
  group_name?: string // 分组名称（如control, treatment）
  layout?: 'PE' | 'SE' | string // 测序类型：双端(PE)或单端(SE)
  batch_id?: string // 批次标识符
  metadata: Record<string, any> // 额外元数据
  created_at: string
  updated_at: string
  file_count?: number // 关联文件数量
}

export interface SampleCreate {
  project_id: string
  sample_id: string
  run_id?: string
  group_name?: string
  layout?: 'PE' | 'SE' | string
  batch_id?: string
  metadata?: Record<string, any>
}

export interface SampleUpdate {
  sample_id?: string
  run_id?: string
  group_name?: string
  layout?: 'PE' | 'SE' | string
  batch_id?: string
  metadata?: Record<string, any>
}

export interface SampleBatchCreate {
  project_id: string
  samples: Array<{
    sample_id: string
    run_id?: string
    group_name?: string
    layout?: 'PE' | 'SE' | string
    batch_id?: string
  }>
}

/**
 * @deprecated Use PaginatedResponse<Sample> instead
 * Kept for backward compatibility during migration
 */
export type SampleListResponse = PaginatedResponse<Sample>
