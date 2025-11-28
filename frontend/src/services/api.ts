/**
 * API Client - Axios wrapper with interceptors
 */
import type { AxiosInstance, AxiosRequestConfig, AxiosResponse } from 'axios'
import axios from 'axios'
import { tokenManager } from '@/utils/tokenManager'

const API_URL = import.meta.env.VITE_API_URL || 'http://localhost:8000/api/v1'

// Flag to prevent multiple redirects
let isRedirecting = false

class ApiClient {
  private client: AxiosInstance

  constructor() {
    this.client = axios.create({
      baseURL: API_URL,
      timeout: 30000,
      headers: {
        'Content-Type': 'application/json',
      },
    })

    this.setupInterceptors()
  }

  private setupInterceptors() {
    // Request interceptor - add auth token
    this.client.interceptors.request.use(
      (config) => {
        const token = tokenManager.getToken()
        if (token) {
          config.headers.Authorization = `Bearer ${token}`
        }
        return config
      },
      (error) => {
        return Promise.reject(error)
      },
    )

    // Response interceptor - handle errors
    this.client.interceptors.response.use(
      (response) => response,
      async (error) => {
        const originalRequest = error.config

        // Handle 401 Unauthorized
        if (error.response?.status === 401 && !originalRequest._retry) {
          originalRequest._retry = true

          // Clear token
          tokenManager.clearToken()

          // Redirect to login (prevent multiple redirects)
          if (!isRedirecting) {
            isRedirecting = true
            // Use a small delay to allow any pending operations to complete
            setTimeout(() => {
              window.location.href = '/login'
            }, 100)
          }

          return Promise.reject(error)
        }

        // Handle other errors
        const message = error.response?.data?.detail || error.message || 'An error occurred'

        return Promise.reject({
          status: error.response?.status,
          message,
          data: error.response?.data,
        })
      },
    )
  }

  async get<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    const response: AxiosResponse<T> = await this.client.get(url, config)
    return response.data
  }

  async post<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const response: AxiosResponse<T> = await this.client.post(url, data, config)
    return response.data
  }

  async put<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const response: AxiosResponse<T> = await this.client.put(url, data, config)
    return response.data
  }

  async delete<T = any>(url: string, config?: AxiosRequestConfig): Promise<T> {
    const response: AxiosResponse<T> = await this.client.delete(url, config)
    return response.data
  }

  async patch<T = any>(url: string, data?: any, config?: AxiosRequestConfig): Promise<T> {
    const response: AxiosResponse<T> = await this.client.patch(url, data, config)
    return response.data
  }
}

export const apiClient = new ApiClient()
export default apiClient
