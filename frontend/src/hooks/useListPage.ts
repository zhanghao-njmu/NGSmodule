/**
 * useListPage Hook
 * Reusable hook for list pages with pagination, search, and filtering
 * Eliminates code duplication across all list page components
 */

import { useState, useEffect, useCallback } from 'react'
import { message } from 'antd'
import type { PaginatedResponse, ListParams } from '@/types/common'

/**
 * Configuration options for useListPage hook
 */
export interface UseListPageOptions<T> {
  /**
   * Function to fetch data from API
   * Should return paginated response
   */
  fetchData: (params?: ListParams) => Promise<PaginatedResponse<T>>

  /**
   * Optional delete function
   */
  deleteItem?: (id: string) => Promise<void>

  /**
   * Initial page number
   * @default 1
   */
  initialPage?: number

  /**
   * Initial page size
   * @default 10
   */
  initialPageSize?: number

  /**
   * Initial search text
   * @default ''
   */
  initialSearch?: string

  /**
   * Initial filters
   * @default {}
   */
  initialFilters?: Record<string, any>

  /**
   * Whether to fetch data automatically on mount
   * @default true
   */
  autoFetch?: boolean

  /**
   * Debounce delay for search (ms)
   * @default 300
   */
  searchDebounce?: number

  /**
   * Error callback
   */
  onError?: (error: Error) => void

  /**
   * Success callback for delete
   */
  onDeleteSuccess?: () => void

  /**
   * Custom error messages
   */
  errorMessages?: {
    fetch?: string
    delete?: string
  }
}

/**
 * Return type for useListPage hook
 */
export interface UseListPageReturn<T> {
  // Data
  items: T[]
  loading: boolean
  total: number
  page: number
  pageSize: number
  searchText: string
  filters: Record<string, any>

  // Actions
  refresh: () => Promise<void>
  handleSearch: (value: string) => void
  handlePageChange: (newPage: number, newPageSize?: number) => void
  handleFilterChange: (filters: Record<string, any>) => void
  handleDelete: (id: string) => Promise<void>
  setItems: React.Dispatch<React.SetStateAction<T[]>>
  setLoading: React.Dispatch<React.SetStateAction<boolean>>
}

/**
 * Hook for list pages with common functionality
 * Handles pagination, search, filtering, and deletion
 *
 * @example
 * const {
 *   items,
 *   loading,
 *   total,
 *   page,
 *   pageSize,
 *   searchText,
 *   refresh,
 *   handleSearch,
 *   handlePageChange,
 *   handleDelete,
 * } = useListPage({
 *   fetchData: projectService.getAll,
 *   deleteItem: projectService.delete,
 * })
 */
export function useListPage<T extends { id: string }>(options: UseListPageOptions<T>): UseListPageReturn<T> {
  const {
    fetchData,
    deleteItem,
    initialPage = 1,
    initialPageSize = 10,
    initialSearch = '',
    initialFilters = {},
    autoFetch = true,
    searchDebounce = 300,
    onError,
    onDeleteSuccess,
    errorMessages = {},
  } = options

  // State
  const [items, setItems] = useState<T[]>([])
  const [loading, setLoading] = useState(false)
  const [total, setTotal] = useState(0)
  const [page, setPage] = useState(initialPage)
  const [pageSize, setPageSize] = useState(initialPageSize)
  const [searchText, setSearchText] = useState(initialSearch)
  const [filters, setFilters] = useState<Record<string, any>>(initialFilters)
  const [searchTimeout, setSearchTimeout] = useState<NodeJS.Timeout | null>(null)

  /**
   * Load data from API
   */
  const loadData = useCallback(
    async (params?: Partial<ListParams>) => {
      setLoading(true)

      try {
        const response = await fetchData({
          page,
          page_size: pageSize,
          search: searchText || undefined,
          ...filters,
          ...params,
        })

        setItems(response.items)
        setTotal(response.total)
      } catch (error: any) {
        const errorMessage = errorMessages.fetch || `Failed to load data: ${error.message}`

        message.error(errorMessage)
        onError?.(error)
      } finally {
        setLoading(false)
      }
    },
    [page, pageSize, searchText, filters, fetchData, errorMessages.fetch, onError],
  )

  /**
   * Refresh data (reload current page)
   */
  const refresh = useCallback(async () => {
    await loadData()
  }, [loadData])

  /**
   * Handle search with debouncing
   */
  const handleSearch = useCallback(
    (value: string) => {
      setSearchText(value)
      setPage(1) // Reset to first page on search

      // Clear existing timeout
      if (searchTimeout) {
        clearTimeout(searchTimeout)
      }

      // Set new timeout for debounced search
      const timeout = setTimeout(() => {
        loadData({ search: value, page: 1 })
      }, searchDebounce)

      setSearchTimeout(timeout)
    },
    [searchTimeout, searchDebounce, loadData],
  )

  /**
   * Handle page change
   */
  const handlePageChange = useCallback(
    (newPage: number, newPageSize?: number) => {
      setPage(newPage)

      if (newPageSize && newPageSize !== pageSize) {
        setPageSize(newPageSize)
        setPage(1) // Reset to first page when page size changes
      }
    },
    [pageSize],
  )

  /**
   * Handle filter change
   */
  const handleFilterChange = useCallback((newFilters: Record<string, any>) => {
    setFilters(newFilters)
    setPage(1) // Reset to first page on filter change
  }, [])

  /**
   * Handle item deletion
   */
  const handleDelete = useCallback(
    async (id: string) => {
      if (!deleteItem) {
        console.warn('Delete function not provided')
        return
      }

      try {
        await deleteItem(id)
        message.success('Deleted successfully')
        onDeleteSuccess?.()

        // Refresh data after deletion
        await loadData()
      } catch (error: any) {
        const errorMessage = errorMessages.delete || `Failed to delete: ${error.message}`

        message.error(errorMessage)
        onError?.(error)
      }
    },
    [deleteItem, loadData, errorMessages.delete, onError, onDeleteSuccess],
  )

  /**
   * Load data on mount and when dependencies change
   */
  useEffect(() => {
    if (autoFetch) {
      loadData()
    }

    // Cleanup timeout on unmount
    return () => {
      if (searchTimeout) {
        clearTimeout(searchTimeout)
      }
    }
  }, [page, pageSize, filters]) // Don't include searchText (handled by debounce)

  return {
    // Data
    items,
    loading,
    total,
    page,
    pageSize,
    searchText,
    filters,

    // Actions
    refresh,
    handleSearch,
    handlePageChange,
    handleFilterChange,
    handleDelete,
    setItems,
    setLoading,
  }
}

export default useListPage
