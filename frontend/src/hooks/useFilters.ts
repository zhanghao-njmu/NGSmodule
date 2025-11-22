/**
 * useFilters Hook - Unified filter state management
 *
 * Eliminates repetitive filter state and change handler logic across list components.
 * Supports multiple filter types: search, select, date range, etc.
 *
 * @example
 * const { filters, setFilter, resetFilters, hasActiveFilters } = useFilters({
 *   search: '',
 *   status: 'all',
 *   dateRange: null,
 * })
 */

import { useState, useCallback, useMemo } from 'react'

export type FilterValue = string | number | boolean | null | undefined | any[]

export interface FiltersState {
  [key: string]: FilterValue
}

export interface UseFiltersOptions<T extends FiltersState> {
  /** Initial filter values */
  initialFilters: T
  /** Callback when filters change */
  onChange?: (filters: T) => void
  /** Filter values that should be considered "empty" (default: '', null, undefined, [], 'all') */
  emptyValues?: FilterValue[]
}

export interface UseFiltersReturn<T extends FiltersState> {
  /** Current filter values */
  filters: T
  /** Set a single filter value */
  setFilter: (key: keyof T, value: FilterValue) => void
  /** Set multiple filter values at once */
  setFilters: (updates: Partial<T>) => void
  /** Reset filters to initial values */
  resetFilters: () => void
  /** Whether any filters are currently active (not at default values) */
  hasActiveFilters: boolean
  /** Get filter value by key */
  getFilter: (key: keyof T) => FilterValue
}

const DEFAULT_EMPTY_VALUES: FilterValue[] = ['', null, undefined, [], 'all']

/**
 * Hook for managing filter state in list components
 */
export function useFilters<T extends FiltersState>(
  options: UseFiltersOptions<T>
): UseFiltersReturn<T> {
  const { initialFilters, onChange, emptyValues = DEFAULT_EMPTY_VALUES } = options

  const [filters, setFiltersState] = useState<T>(initialFilters)

  /**
   * Set a single filter value
   */
  const setFilter = useCallback(
    (key: keyof T, value: FilterValue) => {
      const newFilters = { ...filters, [key]: value }
      setFiltersState(newFilters)
      onChange?.(newFilters)
    },
    [filters, onChange]
  )

  /**
   * Set multiple filter values at once
   */
  const setFilters = useCallback(
    (updates: Partial<T>) => {
      const newFilters = { ...filters, ...updates }
      setFiltersState(newFilters)
      onChange?.(newFilters)
    },
    [filters, onChange]
  )

  /**
   * Reset all filters to initial values
   */
  const resetFilters = useCallback(() => {
    setFiltersState(initialFilters)
    onChange?.(initialFilters)
  }, [initialFilters, onChange])

  /**
   * Get filter value by key
   */
  const getFilter = useCallback(
    (key: keyof T) => {
      return filters[key]
    },
    [filters]
  )

  /**
   * Check if any filters are currently active (not at default/empty values)
   */
  const hasActiveFilters = useMemo(() => {
    return Object.entries(filters).some(([key, value]) => {
      const initialValue = initialFilters[key]

      // If value is different from initial value, it's active
      if (value !== initialValue) {
        // But check if it's an "empty" value
        if (Array.isArray(value)) {
          return value.length > 0
        }
        return !emptyValues.includes(value)
      }

      return false
    })
  }, [filters, initialFilters, emptyValues])

  return {
    filters,
    setFilter,
    setFilters,
    resetFilters,
    hasActiveFilters,
    getFilter,
  }
}

/**
 * Helper function to create a filter change handler
 * Useful for passing to FilterBar or other filter components
 */
export function createFilterChangeHandler<T extends FiltersState>(
  setFilter: (key: keyof T, value: FilterValue) => void
) {
  return (key: string, value: FilterValue) => {
    setFilter(key as keyof T, value)
  }
}

export default useFilters
