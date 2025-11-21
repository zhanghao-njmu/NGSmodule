/**
 * User related types
 */

export interface User {
  id: string
  username: string
  email: string
  full_name: string | null
  role: 'user' | 'admin'
  organization: string | null
  is_active: boolean
  storage_quota: number
  storage_used: number
  created_at: string
}

export interface LoginRequest {
  username: string
  password: string
}

export interface RegisterRequest {
  username: string
  email: string
  password: string
  full_name?: string
  organization?: string
}

export interface TokenResponse {
  access_token: string
  token_type: string
}

export interface UpdateUserRequest {
  full_name?: string
  organization?: string
  email?: string
}
