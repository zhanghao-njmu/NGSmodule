/**
 * User Service - Profile and user management
 */
import apiClient from './api'
import type { User } from '@/types/user'

export interface UserProfile {
  username: string
  email: string
  phone?: string
  avatar?: string
  bio?: string
}

export interface UserStats {
  items: number
  samples: number
  tasks: number
  storageUsed: number
  storageTotal: number
}

export interface UserActivity {
  id: number
  type: 'success' | 'error' | 'warning' | 'info'
  title: string
  description: string
  time: string
  timestamp: string
}

export interface ChangePasswordRequest {
  currentPassword: string
  newPassword: string
}

export interface UserSettings {
  language: string
  timezone: string
  dateFormat: string
  theme: 'light' | 'dark'
}

export interface NotificationSettings {
  emailNotifications: boolean
  pipelineComplete: boolean
  taskFailed: boolean
  systemUpdates: boolean
  weeklyReport: boolean
  browserNotifications: boolean
}

export interface PrivacySettings {
  profileVisible: boolean
  activityVisible: boolean
  twoFactorAuth: boolean
  sessionTimeout: number
}

export interface ApiToken {
  id: string
  name: string
  token: string
  createdAt: string
  lastUsed?: string
  expiresAt?: string
  status: 'active' | 'expired'
}

class UserService {
  /**
   * Get user profile
   */
  async getProfile(): Promise<User> {
    const user = await apiClient.get<User>('/users/me')
    return user
  }

  /**
   * Update user profile
   */
  async updateProfile(data: Partial<UserProfile>): Promise<User> {
    const user = await apiClient.put<User>('/users/me', data)
    return user
  }

  /**
   * Change password
   */
  async changePassword(data: ChangePasswordRequest): Promise<void> {
    await apiClient.post('/users/me/password', data)
  }

  /**
   * Upload avatar
   */
  async uploadAvatar(file: File): Promise<{ avatar_url: string }> {
    const formData = new FormData()
    formData.append('file', file)

    const response = await apiClient.post<{ avatar_url: string }>('/users/me/avatar', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    })

    return response
  }

  /**
   * Get user statistics
   */
  async getUserStats(): Promise<UserStats> {
    const stats = await apiClient.get<UserStats>('/users/me/stats')
    return stats
  }

  /**
   * Get user activity history
   */
  async getUserActivity(limit = 10): Promise<UserActivity[]> {
    const activities = await apiClient.get<UserActivity[]>('/users/me/activity', {
      params: { limit },
    })
    return activities
  }

  /**
   * Delete account
   */
  async deleteAccount(): Promise<void> {
    await apiClient.delete('/users/me')
  }

  /**
   * Get user settings
   */
  async getSettings(): Promise<UserSettings> {
    const settings = await apiClient.get<UserSettings>('/users/me/settings')
    return settings
  }

  /**
   * Update user settings
   */
  async updateSettings(data: Partial<UserSettings>): Promise<UserSettings> {
    const settings = await apiClient.put<UserSettings>('/users/me/settings', data)
    return settings
  }

  /**
   * Get notification settings
   */
  async getNotificationSettings(): Promise<NotificationSettings> {
    const settings = await apiClient.get<NotificationSettings>('/users/me/notifications/settings')
    return settings
  }

  /**
   * Update notification settings
   */
  async updateNotificationSettings(data: Partial<NotificationSettings>): Promise<NotificationSettings> {
    const settings = await apiClient.put<NotificationSettings>('/users/me/notifications/settings', data)
    return settings
  }

  /**
   * Get privacy settings
   */
  async getPrivacySettings(): Promise<PrivacySettings> {
    const settings = await apiClient.get<PrivacySettings>('/users/me/privacy')
    return settings
  }

  /**
   * Update privacy settings
   */
  async updatePrivacySettings(data: Partial<PrivacySettings>): Promise<PrivacySettings> {
    const settings = await apiClient.put<PrivacySettings>('/users/me/privacy', data)
    return settings
  }

  /**
   * Get API tokens
   */
  async getApiTokens(): Promise<ApiToken[]> {
    const tokens = await apiClient.get<ApiToken[]>('/users/me/tokens')
    return tokens
  }

  /**
   * Create API token
   */
  async createApiToken(data: { name: string; description?: string }): Promise<ApiToken> {
    const token = await apiClient.post<ApiToken>('/users/me/tokens', data)
    return token
  }

  /**
   * Delete API token
   */
  async deleteApiToken(tokenId: string): Promise<void> {
    await apiClient.delete(`/users/me/tokens/${tokenId}`)
  }
}

export const userService = new UserService()
export default userService
