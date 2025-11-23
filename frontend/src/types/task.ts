/**
 * Task types and interfaces
 */

import type { PaginatedResponse } from './common'

export type TaskStatus = 'pending' | 'running' | 'completed' | 'failed' | 'cancelled'

export interface Task {
  id: string
  project_id: string
  task_name: string
  task_type?: string
  status: TaskStatus
  progress: number
  started_at?: string
  completed_at?: string
  error_message?: string
  config: Record<string, any>
  celery_task_id?: string
  log_file_path?: string
  created_at: string
}

export interface TaskCreate {
  project_id: string
  task_name: string
  task_type?: string
  config?: Record<string, any>
}

export interface TaskUpdate {
  task_name?: string
  task_type?: string
  status?: TaskStatus
  progress?: number
  error_message?: string
  config?: Record<string, any>
}

export interface TaskExecuteRequest {
  pipeline_script: string
  config?: Record<string, any>
}

/**
 * @deprecated Use PaginatedResponse<Task> instead
 */
export type TaskListResponse = PaginatedResponse<Task>

export interface TaskStats {
  total_tasks: number
  pending_tasks: number
  running_tasks: number
  completed_tasks: number
  failed_tasks: number
  cancelled_tasks: number
}

export interface TaskLogResponse {
  task_id: string
  log_content: string
  log_file_path?: string
}

export interface TaskUpdateEvent {
  type: 'task_update'
  task_id: string
  status: TaskStatus
  progress: number
  message: string
  timestamp: string
}

export interface WebSocketMessage {
  type: 'subscribe' | 'unsubscribe' | 'ping' | 'task_update' | 'subscribed' | 'unsubscribed' | 'pong' | 'error'
  task_id?: string
  status?: TaskStatus
  progress?: number
  message?: string
  timestamp?: string
}
