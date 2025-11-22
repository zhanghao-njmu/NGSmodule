/**
 * CRUD Service Factory
 * Generic factory for creating CRUD (Create, Read, Update, Delete) services
 * Eliminates code duplication across service layers
 */

import apiClient from './api'
import { PaginatedResponse, ListParams } from '@/types/common'

/**
 * Generic CRUD service interface
 * Defines standard CRUD operations for any entity
 *
 * @template T - The entity type
 * @template CreateT - The create payload type (defaults to Partial<T>)
 * @template UpdateT - The update payload type (defaults to Partial<T>)
 */
export interface CrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>> {
  /**
   * Get all items with optional filters/pagination
   */
  getAll: (params?: ListParams) => Promise<PaginatedResponse<T>>

  /**
   * Get a single item by ID
   */
  getById: (id: string) => Promise<T>

  /**
   * Create a new item
   */
  create: (data: CreateT) => Promise<T>

  /**
   * Update an existing item
   */
  update: (id: string, data: UpdateT) => Promise<T>

  /**
   * Delete an item by ID
   */
  delete: (id: string) => Promise<void>
}

/**
 * Extended CRUD service interface with additional common methods
 */
export interface ExtendedCrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>>
  extends CrudService<T, CreateT, UpdateT> {
  /**
   * Batch delete multiple items
   */
  batchDelete?: (ids: string[]) => Promise<void>

  /**
   * Duplicate/clone an existing item
   */
  duplicate?: (id: string) => Promise<T>

  /**
   * Export items to file
   */
  export?: (params?: ListParams) => Promise<{ downloadUrl: string }>

  /**
   * Import items from file
   */
  import?: (file: File) => Promise<{ imported: number; failed: number }>
}

/**
 * Configuration options for CRUD service factory
 */
export interface CrudServiceConfig {
  /**
   * API endpoint (e.g., 'projects', 'samples')
   * Will be prefixed with '/' automatically
   */
  endpoint: string

  /**
   * Whether the API returns paginated responses
   * If false, will wrap non-paginated responses in PaginatedResponse format
   * @default true
   */
  paginated?: boolean

  /**
   * Transform function for list responses
   * Useful for adapting backend response format
   */
  transformList?: (response: any) => PaginatedResponse<any>

  /**
   * Transform function for single item responses
   */
  transformItem?: (response: any) => any

  /**
   * Custom headers for requests
   */
  headers?: Record<string, string>
}

/**
 * Create a CRUD service for a given entity type
 *
 * @param config - Service configuration
 * @returns CRUD service instance
 *
 * @example
 * // Basic usage
 * const projectService = createCrudService<Project>({ endpoint: 'projects' })
 *
 * // With custom types
 * const sampleService = createCrudService<Sample, SampleCreate, SampleUpdate>({
 *   endpoint: 'samples'
 * })
 *
 * // With transformation
 * const taskService = createCrudService<Task>({
 *   endpoint: 'tasks',
 *   transformList: (response) => ({
 *     total: response.count,
 *     items: response.results,
 *   })
 * })
 */
export function createCrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>>(
  config: CrudServiceConfig
): CrudService<T, CreateT, UpdateT> {
  const { endpoint, paginated = true, transformList, transformItem, headers } = config

  // Ensure endpoint starts with /
  const baseEndpoint = endpoint.startsWith('/') ? endpoint : `/${endpoint}`

  return {
    /**
     * Get all items with optional filters and pagination
     */
    getAll: async (params?: ListParams): Promise<PaginatedResponse<T>> => {
      const response = await apiClient.get<any>(baseEndpoint, {
        params,
        headers,
      })

      // Transform response if transformer provided
      if (transformList) {
        return transformList(response)
      }

      // If response is already in paginated format
      if (paginated && response.items && typeof response.total === 'number') {
        return response as PaginatedResponse<T>
      }

      // Wrap non-paginated response
      return {
        total: Array.isArray(response) ? response.length : response.total || 0,
        items: Array.isArray(response) ? response : response.items || [],
      }
    },

    /**
     * Get a single item by ID
     */
    getById: async (id: string): Promise<T> => {
      const response = await apiClient.get<T>(`${baseEndpoint}/${id}`, { headers })

      if (transformItem) {
        return transformItem(response)
      }

      return response
    },

    /**
     * Create a new item
     */
    create: async (data: CreateT): Promise<T> => {
      const response = await apiClient.post<T>(baseEndpoint, data, { headers })

      if (transformItem) {
        return transformItem(response)
      }

      return response
    },

    /**
     * Update an existing item
     */
    update: async (id: string, data: UpdateT): Promise<T> => {
      const response = await apiClient.put<T>(`${baseEndpoint}/${id}`, data, { headers })

      if (transformItem) {
        return transformItem(response)
      }

      return response
    },

    /**
     * Delete an item by ID
     */
    delete: async (id: string): Promise<void> => {
      await apiClient.delete(`${baseEndpoint}/${id}`, { headers })
    },
  }
}

/**
 * Create an extended CRUD service with additional common operations
 *
 * @param config - Service configuration
 * @returns Extended CRUD service instance
 *
 * @example
 * const projectService = createExtendedCrudService<Project>({
 *   endpoint: 'projects'
 * })
 *
 * // Now supports batch delete, duplicate, export, import
 * await projectService.batchDelete(['id1', 'id2'])
 * await projectService.duplicate('project-id')
 * const exportUrl = await projectService.export({ status: 'completed' })
 */
export function createExtendedCrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>>(
  config: CrudServiceConfig
): ExtendedCrudService<T, CreateT, UpdateT> {
  const baseService = createCrudService<T, CreateT, UpdateT>(config)
  const { endpoint, headers } = config
  const baseEndpoint = endpoint.startsWith('/') ? endpoint : `/${endpoint}`

  return {
    ...baseService,

    /**
     * Batch delete multiple items
     */
    batchDelete: async (ids: string[]): Promise<void> => {
      await apiClient.post(
        `${baseEndpoint}/batch-delete`,
        { ids },
        { headers }
      )
    },

    /**
     * Duplicate an existing item
     */
    duplicate: async (id: string): Promise<T> => {
      const response = await apiClient.post<T>(
        `${baseEndpoint}/${id}/duplicate`,
        {},
        { headers }
      )

      if (config.transformItem) {
        return config.transformItem(response)
      }

      return response
    },

    /**
     * Export items to downloadable file
     */
    export: async (params?: ListParams): Promise<{ downloadUrl: string }> => {
      const response = await apiClient.get<{ download_url: string }>(
        `${baseEndpoint}/export`,
        { params, headers }
      )

      return {
        downloadUrl: response.download_url,
      }
    },

    /**
     * Import items from file
     */
    import: async (file: File): Promise<{ imported: number; failed: number }> => {
      const formData = new FormData()
      formData.append('file', file)

      const response = await apiClient.post<{ imported: number; failed: number }>(
        `${baseEndpoint}/import`,
        formData,
        {
          headers: {
            'Content-Type': 'multipart/form-data',
            ...headers,
          },
        }
      )

      return response
    },
  }
}

/**
 * Utility: Extend a CRUD service with custom methods
 *
 * @param baseService - Base CRUD service
 * @param customMethods - Additional methods to add
 * @returns Extended service with custom methods
 *
 * @example
 * const baseService = createCrudService<Project>({ endpoint: 'projects' })
 *
 * const projectService = extendService(baseService, {
 *   async getStats(projectId: string) {
 *     return apiClient.get(`/projects/${projectId}/stats`)
 *   },
 *   async archive(projectId: string) {
 *     return apiClient.post(`/projects/${projectId}/archive`)
 *   }
 * })
 */
export function extendService<T extends Record<string, any>, E extends Record<string, any>>(
  baseService: T,
  customMethods: E
): T & E {
  return {
    ...baseService,
    ...customMethods,
  }
}

/**
 * Batch create multiple items
 * Helper function for services that need batch creation
 *
 * @param endpoint - API endpoint
 * @param items - Array of items to create
 * @returns Array of created items
 *
 * @example
 * await batchCreate('/samples', [sample1, sample2, sample3])
 */
export async function batchCreate<T, CreateT = Partial<T>>(
  endpoint: string,
  items: CreateT[]
): Promise<T[]> {
  const response = await apiClient.post<T[]>(
    `${endpoint}/batch`,
    { items }
  )

  return response
}

/**
 * Batch update multiple items
 * Helper function for services that need batch updates
 *
 * @param endpoint - API endpoint
 * @param updates - Array of updates with IDs
 * @returns Array of updated items
 *
 * @example
 * await batchUpdate('/samples', [
 *   { id: '1', data: { status: 'active' } },
 *   { id: '2', data: { status: 'inactive' } }
 * ])
 */
export async function batchUpdate<T, UpdateT = Partial<T>>(
  endpoint: string,
  updates: Array<{ id: string; data: UpdateT }>
): Promise<T[]> {
  const response = await apiClient.put<T[]>(
    `${endpoint}/batch`,
    { updates }
  )

  return response
}

/**
 * Search items with full-text search
 * Helper function for services with search functionality
 *
 * @param endpoint - API endpoint
 * @param query - Search query
 * @param params - Additional search parameters
 * @returns Paginated search results
 *
 * @example
 * const results = await searchItems('/projects', 'RNA-seq', { status: 'active' })
 */
export async function searchItems<T>(
  endpoint: string,
  query: string,
  params?: Record<string, any>
): Promise<PaginatedResponse<T>> {
  const response = await apiClient.get<PaginatedResponse<T>>(
    `${endpoint}/search`,
    {
      params: {
        q: query,
        ...params,
      },
    }
  )

  return response
}

export default createCrudService
