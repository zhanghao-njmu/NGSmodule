/**
 * usePagination Hook - Unified pagination state management
 */
import { useState, useCallback } from 'react'

export interface PaginationState {
  current: number
  pageSize: number
  total: number
}

export interface UsePaginationOptions {
  initialPage?: number
  initialPageSize?: number
  initialTotal?: number
}

export function usePagination(options: UsePaginationOptions = {}) {
  const {
    initialPage = 1,
    initialPageSize = 20,
    initialTotal = 0,
  } = options

  const [pagination, setPagination] = useState<PaginationState>({
    current: initialPage,
    pageSize: initialPageSize,
    total: initialTotal,
  })

  const setPage = useCallback((page: number) => {
    setPagination((prev) => ({ ...prev, current: page }))
  }, [])

  const setPageSize = useCallback((size: number) => {
    setPagination((prev) => ({ ...prev, pageSize: size, current: 1 }))
  }, [])

  const setTotal = useCallback((total: number) => {
    setPagination((prev) => ({ ...prev, total }))
  }, [])

  const reset = useCallback(() => {
    setPagination({
      current: initialPage,
      pageSize: initialPageSize,
      total: initialTotal,
    })
  }, [initialPage, initialPageSize, initialTotal])

  const onChange = useCallback((page: number, pageSize?: number) => {
    setPagination((prev) => ({
      ...prev,
      current: page,
      pageSize: pageSize || prev.pageSize,
    }))
  }, [])

  // Calculate skip for API calls
  const skip = (pagination.current - 1) * pagination.pageSize

  return {
    pagination,
    setPage,
    setPageSize,
    setTotal,
    reset,
    onChange,
    skip,
    limit: pagination.pageSize,
  }
}

export default usePagination
