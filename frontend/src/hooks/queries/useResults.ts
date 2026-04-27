/**
 * TanStack Query hooks for the Results domain.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import resultService from '@/services/result.service'

export function useResultList(params?: any) {
  return useQuery({
    queryKey: [...queryKeys.results.all, 'list', params],
    queryFn: () => resultService.getAll(params),
  })
}

export function useResultsByTask(taskId: string | undefined, params?: any) {
  return useQuery({
    queryKey: queryKeys.results.list(taskId ?? ''),
    queryFn: () => resultService.getResultsByTask(taskId as string, params),
    enabled: !!taskId,
  })
}

export function useResultsByType(type: string | undefined, params?: any) {
  return useQuery({
    queryKey: [...queryKeys.results.all, 'by-type', type, params],
    queryFn: () => resultService.getResultsByType(type as string, params),
    enabled: !!type,
  })
}

export function useResult(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.results.detail(id ?? ''),
    queryFn: () => resultService.getById(id as string),
    enabled: !!id,
  })
}

export function useResultVisualization(id: string | undefined) {
  return useQuery({
    queryKey: [...queryKeys.results.detail(id ?? ''), 'viz'],
    queryFn: () => resultService.getVisualizationData(id as string),
    enabled: !!id,
    // Visualization data is derived/expensive — cache longer
    staleTime: 5 * 60_000,
  })
}

export function useTaskResultsSummary(taskId: string | undefined) {
  return useQuery({
    queryKey: [...queryKeys.results.all, 'task-summary', taskId],
    queryFn: () => resultService.getTaskResultsSummary(taskId as string),
    enabled: !!taskId,
  })
}

// ----- mutations -----

export function useExportResult() {
  return useMutation({
    mutationFn: ({ id, format }: { id: string; format?: 'csv' | 'json' | 'tsv' }) =>
      resultService.exportResult(id, format),
  })
}

export function useDeleteResult() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => resultService.delete(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.results.all })
    },
  })
}
