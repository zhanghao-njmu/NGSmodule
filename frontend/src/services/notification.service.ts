/**
 * Notification Service
 */
import apiClient from './api'
import type {
  Notification,
  NotificationStats,
  NotificationListParams,
} from '@/types/notification'

export interface NotificationListResponse {
  notifications: Notification[]
  total: number
  page: number
  pageSize: number
}

class NotificationService {
  /**
   * Get notifications list
   */
  async getNotifications(
    params?: NotificationListParams
  ): Promise<NotificationListResponse> {
    const response = await apiClient.get<NotificationListResponse>('/notifications', {
      params,
    })
    return response
  }

  /**
   * Get notification by ID
   */
  async getNotification(id: string): Promise<Notification> {
    const notification = await apiClient.get<Notification>(`/notifications/${id}`)
    return notification
  }

  /**
   * Mark notification as read
   */
  async markAsRead(id: string): Promise<void> {
    await apiClient.put(`/notifications/${id}/read`)
  }

  /**
   * Mark multiple notifications as read
   */
  async markMultipleAsRead(ids: string[]): Promise<void> {
    await apiClient.put('/notifications/read', { ids })
  }

  /**
   * Mark all notifications as read
   */
  async markAllAsRead(): Promise<void> {
    await apiClient.put('/notifications/read-all')
  }

  /**
   * Delete notification
   */
  async deleteNotification(id: string): Promise<void> {
    await apiClient.delete(`/notifications/${id}`)
  }

  /**
   * Delete multiple notifications
   */
  async deleteMultiple(ids: string[]): Promise<void> {
    await apiClient.delete('/notifications', { data: { ids } })
  }

  /**
   * Delete all read notifications
   */
  async deleteAllRead(): Promise<void> {
    await apiClient.delete('/notifications/read')
  }

  /**
   * Get notification statistics
   */
  async getStats(): Promise<NotificationStats> {
    const stats = await apiClient.get<NotificationStats>('/notifications/stats')
    return stats
  }

  /**
   * Get unread count
   */
  async getUnreadCount(): Promise<number> {
    const { count } = await apiClient.get<{ count: number }>('/notifications/unread-count')
    return count
  }
}

export const notificationService = new NotificationService()
export default notificationService
