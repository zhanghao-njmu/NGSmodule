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

    const response = await apiClient.post<{ avatar_url: string }>(
      '/users/me/avatar',
      formData,
      {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      }
    )

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
}

export const userService = new UserService()
export default userService
