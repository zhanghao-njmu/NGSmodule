/**
 * Project types and interfaces
 */

import type { PaginatedResponse } from './common'

export interface Project {
  id: string
  name: string
  description?: string
  user_id: string
  status: 'active' | 'archived' | 'completed'
  config: Record<string, any>
  created_at: string
  updated_at: string
  sample_count?: number
  task_count?: number
}

export interface ProjectCreate {
  name: string
  description?: string
  config?: Record<string, any>
}

export interface ProjectUpdate {
  name?: string
  description?: string
  status?: 'active' | 'archived' | 'completed'
  config?: Record<string, any>
}

/**
 * @deprecated Use PaginatedResponse<Project> instead
 * Kept for backward compatibility during migration
 */
export type ProjectListResponse = PaginatedResponse<Project>

export interface ProjectStats {
  total_projects: number
  active_projects: number
  archived_projects: number
  completed_projects: number
  total_tasks: number
  active_tasks: number
}
