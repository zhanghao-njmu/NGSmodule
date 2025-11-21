/**
 * Admin types and interfaces
 */

export interface UserStats {
  user_id: string
  username: string
  total_projects: number
  total_samples: number
  total_tasks: number
  completed_tasks: number
  failed_tasks: number
  storage_used: number
  storage_quota: number
  storage_percent: number
}

export interface SystemStats {
  total_users: number
  active_users: number
  total_projects: number
  total_samples: number
  total_tasks: number
  running_tasks: number
  completed_tasks: number
  failed_tasks: number
  total_storage_used: number
  total_storage_quota: number
}

export interface User {
  id: string
  username: string
  email: string
  full_name?: string
  role: string
  organization?: string
  is_active: boolean
  storage_quota: number
  storage_used: number
  created_at: string
}

export interface UserAdminUpdate {
  full_name?: string
  organization?: string
  email?: string
  role?: 'user' | 'admin'
  is_active?: boolean
  storage_quota?: number
}
