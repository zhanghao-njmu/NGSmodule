/**
 * Notification Types
 */

export type NotificationType = 'success' | 'info' | 'warning' | 'error'

export type NotificationCategory = 'pipeline' | 'task' | 'system' | 'security' | 'project' | 'sample' | 'result'

export interface Notification {
  id: string
  type: NotificationType
  category: NotificationCategory
  title: string
  message: string
  timestamp: string
  read: boolean
  actionUrl?: string
  actionText?: string
  metadata?: Record<string, any>
}

export interface NotificationStats {
  total: number
  unread: number
  byCategory: Record<NotificationCategory, number>
  byType: Record<NotificationType, number>
}

export interface NotificationFilter {
  type?: NotificationType[]
  category?: NotificationCategory[]
  read?: boolean
  startDate?: string
  endDate?: string
}

export interface NotificationListParams {
  page?: number
  pageSize?: number
  filter?: NotificationFilter
  sortBy?: 'timestamp' | 'type' | 'category'
  sortOrder?: 'asc' | 'desc'
}
