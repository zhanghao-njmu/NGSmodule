/**
 * TanStack Query hooks for the Tasks domain.
 *
 * Adapted to the existing task service which uses the CRUD factory.
 *
 * Polling strategy:
 *   - List: 10s stale window, refetch on window focus
 *   - Detail: 5s polling while running, stop once terminal
 *   - Realtime layer (useTaskProgress) is the primary live update path;
 *     polling acts as a safety net.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import taskService from '@/services/task.service'

export interface TaskListParams {
  project_id?: string
  status?: string
  task_type?: string
  skip?: number
  limit?: number
}

const TERMINAL_STATUSES = new Set(['completed', 'failed', 'cancelled'])

export function useTaskList(params: TaskListParams = {}) {
  return useQuery({
    queryKey: queryKeys.tasks.list(params),
    queryFn: () => taskService.getAll(params as any),
    staleTime: 10_000,
  })
}

export function useTaskDetail(
  taskId: string | undefined,
  options: { enablePolling?: boolean } = { enablePolling: true },
) {
  return useQuery({
    queryKey: queryKeys.tasks.detail(taskId ?? ''),
    queryFn: () => taskService.getById(taskId as string),
    enabled: !!taskId,
    refetchInterval: (query) => {
      if (!options.enablePolling) return false
      const status = (query.state.data as any)?.status
      if (status && TERMINAL_STATUSES.has(status)) return false
      return 5_000
    },
  })
}

export function useTaskLogs(taskId: string | undefined) {
  return useQuery({
    queryKey: queryKeys.tasks.logs(taskId ?? ''),
    queryFn: () => taskService.getTaskLogs(taskId as string),
    enabled: !!taskId,
    staleTime: 5_000,
  })
}

export function useTaskStatsSummary(params?: { project_id?: string }) {
  return useQuery({
    queryKey: [...queryKeys.tasks.stats, params],
    queryFn: () => taskService.getStats(params),
    staleTime: 30_000,
  })
}

// ----- mutations -----

export function useCreateTask() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => taskService.create(data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.stats.tasks })
    },
  })
}

export function useExecuteTask() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ taskId, data }: { taskId: string; data: any }) =>
      taskService.executeTask(taskId, data),
    onSuccess: (_, { taskId }) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.detail(taskId) })
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
    },
  })
}

export function useCancelTask() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (taskId: string) => taskService.cancelTask(taskId),
    onSuccess: (_, taskId) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.detail(taskId) })
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
    },
  })
}

export function useRetryTask() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (taskId: string) => taskService.retryTask(taskId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
    },
  })
}

export function useDeleteTask() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (taskId: string) => taskService.delete(taskId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })
    },
  })
}
