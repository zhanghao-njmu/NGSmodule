/**
 * useAsync Hook - Unified async operation handling
 *
 * Eliminates repetitive loading/error/data state management
 */
import { useState, useCallback, useEffect } from 'react'

export interface AsyncState<T> {
  data: T | null
  loading: boolean
  error: Error | null
}

export interface UseAsyncOptions {
  immediate?: boolean // Execute on mount
  onSuccess?: (data: any) => void
  onError?: (error: Error) => void
}

export function useAsync<T>(asyncFunction: (...args: any[]) => Promise<T>, options: UseAsyncOptions = {}) {
  const { immediate = false, onSuccess, onError } = options

  const [state, setState] = useState<AsyncState<T>>({
    data: null,
    loading: false,
    error: null,
  })

  const execute = useCallback(
    async (...args: any[]) => {
      setState({ data: null, loading: true, error: null })

      try {
        const data = await asyncFunction(...args)
        setState({ data, loading: false, error: null })
        onSuccess?.(data)
        return data
      } catch (error) {
        const err = error instanceof Error ? error : new Error(String(error))
        setState({ data: null, loading: false, error: err })
        onError?.(err)
        throw err
      }
    },
    [asyncFunction, onSuccess, onError],
  )

  // Execute immediately on mount if requested
  useEffect(() => {
    if (immediate) {
      execute()
    }
  }, [immediate, execute])

  const reset = useCallback(() => {
    setState({ data: null, loading: false, error: null })
  }, [])

  return {
    ...state,
    execute,
    reset,
  }
}

export default useAsync
