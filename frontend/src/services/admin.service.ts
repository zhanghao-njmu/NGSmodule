/**
 * Admin Service - API calls for admin management
 */
import apiClient from './api'
import type { User, UserAdminUpdate, UserStats, SystemStats } from '../types/admin'

export const adminService = {
  /**
   * Get all users
   */
  async getUsers(params?: { skip?: number; limit?: number }): Promise<User[]> {
    return apiClient.get('/users', { params })
  },

  /**
   * Get user by ID
   */
  async getUser(userId: string): Promise<User> {
    return apiClient.get(`/users/${userId}`)
  },

  /**
   * Update user (admin)
   */
  async updateUser(userId: string, data: UserAdminUpdate): Promise<User> {
    return apiClient.put(`/users/${userId}`, data)
  },

  /**
   * Toggle user active status
   */
  async toggleUserStatus(userId: string): Promise<{ message: string }> {
    return apiClient.post(`/users/${userId}/toggle`)
  },

  /**
   * Delete user
   */
  async deleteUser(userId: string): Promise<void> {
    return apiClient.delete(`/users/${userId}`)
  },

  /**
   * Get user statistics
   */
  async getUserStats(userId: string): Promise<UserStats> {
    return apiClient.get(`/users/${userId}/stats`)
  },

  /**
   * Get system statistics
   */
  async getSystemStats(): Promise<SystemStats> {
    return apiClient.get('/users/stats/system')
  },
}

export default adminService
