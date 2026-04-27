/**
 * Notification Service
 * 通知服务 - 对接后端 /api/v1/notifications/* 端点
 */
import apiClient from './api'
import type { Notification, NotificationStats, NotificationListParams } from '@/types/notification'

export interface NotificationListResponse {
  items: Notification[]
  total: number
  page: number
  page_size: number
  unread_count: number
}

export interface UnreadCount {
  count: number
  by_type: Record<string, number>
}

export interface NotificationSettings {
  id?: string
  user_id?: string
  email_enabled: boolean
  email_task_completed?: boolean
  email_task_failed?: boolean
  email_system_alerts?: boolean
  app_enabled: boolean
  app_task_updates?: boolean
  app_project_updates?: boolean
  app_system_alerts?: boolean
  push_enabled: boolean
  updated_at?: string
}

export interface MarkAllReadResponse {
  marked_count: number
  message: string
}

class NotificationService {
  /**
   * 获取通知列表（分页）
   */
  async getNotifications(params?: NotificationListParams): Promise<NotificationListResponse> {
    return apiClient.get<NotificationListResponse>('/notifications', { params })
  }

  /**
   * 获取通知详情
   */
  async getNotification(id: string): Promise<Notification> {
    return apiClient.get<Notification>(`/notifications/${id}`)
  }

  /**
   * 标记通知为已读
   */
  async markAsRead(id: string): Promise<Notification> {
    return apiClient.put<Notification>(`/notifications/${id}/read`)
  }

  /**
   * 标记多个通知为已读
   * 注意: 后端目前不支持批量标记，需逐个调用
   */
  async markMultipleAsRead(ids: string[]): Promise<void> {
    await Promise.all(ids.map((id) => this.markAsRead(id)))
  }

  /**
   * 标记所有通知为已读
   */
  async markAllAsRead(): Promise<MarkAllReadResponse> {
    return apiClient.put<MarkAllReadResponse>('/notifications/read-all')
  }

  /**
   * 删除单个通知
   */
  async deleteNotification(id: string): Promise<void> {
    await apiClient.delete(`/notifications/${id}`)
  }

  /**
   * 删除多个通知（前端批量调用）
   */
  async deleteMultiple(ids: string[]): Promise<void> {
    await Promise.all(ids.map((id) => this.deleteNotification(id)))
  }

  /**
   * 删除所有已读通知（前端实现）
   */
  async deleteAllRead(): Promise<void> {
    const list = await this.getNotifications({ page: 1, page_size: 1000 } as NotificationListParams)
    const readIds = list.items.filter((n: any) => n.read).map((n: any) => n.id)
    if (readIds.length > 0) {
      await this.deleteMultiple(readIds)
    }
  }

  /**
   * 获取通知统计（兼容旧接口，使用未读计数）
   */
  async getStats(): Promise<NotificationStats> {
    const unread = await this.getUnreadCount()
    return {
      total: 0, // 需要单独获取
      unread: unread.count,
      by_type: unread.by_type,
    } as unknown as NotificationStats
  }

  /**
   * 获取未读通知数量（含按类型分组）
   */
  async getUnreadCount(): Promise<UnreadCount> {
    return apiClient.get<UnreadCount>('/notifications/unread/count')
  }

  /**
   * 获取仅未读数字（兼容旧接口）
   */
  async getUnreadCountNumber(): Promise<number> {
    const result = await this.getUnreadCount()
    return result.count
  }

  /**
   * 获取当前用户的通知设置
   */
  async getSettings(): Promise<NotificationSettings> {
    return apiClient.get<NotificationSettings>('/notifications/settings/current')
  }

  /**
   * 更新当前用户的通知设置
   */
  async updateSettings(settings: Partial<NotificationSettings>): Promise<NotificationSettings> {
    return apiClient.put<NotificationSettings>('/notifications/settings/current', settings)
  }
}

export const notificationService = new NotificationService()
export default notificationService
