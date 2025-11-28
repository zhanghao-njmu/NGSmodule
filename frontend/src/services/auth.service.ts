/**
 * Authentication Service
 */
import apiClient from './api'
import { tokenManager } from '@/utils/tokenManager'
import type { LoginRequest, RegisterRequest, TokenResponse, User } from '@/types/user'

class AuthService {
  /**
   * Login user
   */
  async login(credentials: LoginRequest): Promise<{ token: string; user: User }> {
    // FastAPI OAuth2 expects form data
    const formData = new FormData()
    formData.append('username', credentials.username)
    formData.append('password', credentials.password)

    const tokenResponse = await apiClient.post<TokenResponse>('/auth/login', formData, {
      headers: {
        'Content-Type': 'multipart/form-data',
      },
    })

    const token = tokenResponse.access_token

    // Store token using secure token manager
    tokenManager.setToken(token)

    // Get user info
    const user = await this.getCurrentUser()

    return { token, user }
  }

  /**
   * Register new user
   */
  async register(userData: RegisterRequest): Promise<User> {
    const user = await apiClient.post<User>('/auth/register', userData)
    return user
  }

  /**
   * Logout user
   */
  async logout(): Promise<void> {
    try {
      await apiClient.post('/auth/logout')
    } finally {
      tokenManager.clearToken()
    }
  }

  /**
   * Get current user info
   */
  async getCurrentUser(): Promise<User> {
    const user = await apiClient.get<User>('/users/me')
    return user
  }

  /**
   * Check if user is authenticated
   */
  isAuthenticated(): boolean {
    return tokenManager.hasValidToken()
  }

  /**
   * Get stored token
   */
  getToken(): string | null {
    return tokenManager.getToken()
  }
}

export const authService = new AuthService()
export default authService
