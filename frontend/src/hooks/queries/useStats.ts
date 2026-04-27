/**
 * TanStack Query hooks for the Stats / Dashboard domain.
 */
import { useQuery } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import statsService from '@/services/stats.service'

export function useQuickStats() {
  return useQuery({
    queryKey: queryKeys.stats.quick,
    queryFn: () => statsService.getQuickStats(),
    staleTime: 30_000,
    refetchInterval: 60_000,
  })
}

export function useStatsSummary() {
  return useQuery({
    queryKey: queryKeys.stats.summary,
    queryFn: () => statsService.getSummary(),
    staleTime: 30_000,
  })
}

export function useProjectStats() {
  return useQuery({
    queryKey: queryKeys.stats.projects,
    queryFn: () => statsService.getProjectStats(),
  })
}

export function useSampleStats() {
  return useQuery({
    queryKey: queryKeys.stats.samples,
    queryFn: () => statsService.getSampleStats(),
  })
}

export function useTaskStats() {
  return useQuery({
    queryKey: queryKeys.stats.tasks,
    queryFn: () => statsService.getTaskStats(),
  })
}

export function useStorageStats() {
  return useQuery({
    queryKey: queryKeys.stats.storage,
    queryFn: () => statsService.getStorageStats(),
  })
}

export function useTrends(
  metric: 'tasks' | 'samples' | 'storage' | 'projects',
  period: 'hourly' | 'daily' | 'weekly' | 'monthly' = 'daily',
  days = 30,
) {
  return useQuery({
    queryKey: queryKeys.stats.trends(metric, period),
    queryFn: () => statsService.getTrends(metric, period, days),
    staleTime: 5 * 60_000,
  })
}

/**
 * Dashboard composite hook — wraps the legacy `getDashboardStats` shim
 * for compatibility, but everything underneath benefits from query
 * caching and realtime invalidation.
 */
export function useDashboardStats() {
  return useQuery({
    queryKey: ['stats', 'dashboard'],
    queryFn: () => statsService.getDashboardStats(),
    staleTime: 30_000,
    refetchInterval: 60_000,
  })
}
