/**
 * TanStack Query client configuration.
 *
 * Tuned for NGSmodule's mix of:
 *   - Long-running pipeline tasks (frequent refetch)
 *   - Mostly-static reference data (templates, knowledge base)
 *   - Real-time notifications (cache invalidated by WebSocket events)
 */
import { QueryClient } from '@tanstack/react-query'
import { message } from 'antd'

export const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      // Default: keep server state fresh for 30s, then revalidate in background
      staleTime: 30_000,
      // Keep cached data for 5 minutes after last component unmount
      gcTime: 5 * 60_000,
      retry: (failureCount, error: any) => {
        // Don't retry on 4xx (client errors)
        const status = error?.response?.status
        if (status && status >= 400 && status < 500) {
          return false
        }
        return failureCount < 2
      },
      // Refetch when window regains focus (catches background changes)
      refetchOnWindowFocus: true,
      // Don't refetch on mount if we have fresh data
      refetchOnMount: 'always',
    },
    mutations: {
      retry: false,
      onError: (error: any) => {
        const detail = error?.response?.data?.detail || error?.message || 'Operation failed'
        // Surface mutation errors as notifications by default; individual
        // mutations can override with their own onError
        if (typeof detail === 'string') {
          message.error(detail)
        }
      },
    },
  },
})

/**
 * Centralized cache key factories.
 * Using factories prevents typos and makes invalidation predictable.
 */
export const queryKeys = {
  // Auth / user
  auth: {
    me: ['auth', 'me'] as const,
  },
  // Projects
  projects: {
    all: ['projects'] as const,
    list: (params?: Record<string, any>) => ['projects', 'list', params] as const,
    detail: (id: string) => ['projects', 'detail', id] as const,
    stats: (id: string) => ['projects', 'stats', id] as const,
  },
  // Samples
  samples: {
    all: ['samples'] as const,
    list: (projectId: string, params?: Record<string, any>) => ['samples', 'list', projectId, params] as const,
    detail: (id: string) => ['samples', 'detail', id] as const,
  },
  // Tasks
  tasks: {
    all: ['tasks'] as const,
    list: (params?: Record<string, any>) => ['tasks', 'list', params] as const,
    detail: (id: string) => ['tasks', 'detail', id] as const,
    logs: (id: string) => ['tasks', 'logs', id] as const,
    stats: ['tasks', 'stats'] as const,
  },
  // Files
  files: {
    all: ['files'] as const,
    list: (params?: Record<string, any>) => ['files', 'list', params] as const,
    detail: (id: string) => ['files', 'detail', id] as const,
  },
  // Notifications
  notifications: {
    all: ['notifications'] as const,
    list: (params?: Record<string, any>) => ['notifications', 'list', params] as const,
    unreadCount: ['notifications', 'unread-count'] as const,
    settings: ['notifications', 'settings'] as const,
  },
  // Pipelines
  pipelines: {
    all: ['pipelines'] as const,
    templates: ['pipelines', 'templates'] as const,
    template: (id: string) => ['pipelines', 'template', id] as const,
  },
  // Results
  results: {
    all: ['results'] as const,
    list: (taskId: string) => ['results', 'list', taskId] as const,
    detail: (id: string) => ['results', 'detail', id] as const,
  },
  // Stats / analytics
  stats: {
    all: ['stats'] as const,
    summary: ['stats', 'summary'] as const,
    quick: ['stats', 'quick'] as const,
    projects: ['stats', 'projects'] as const,
    samples: ['stats', 'samples'] as const,
    tasks: ['stats', 'tasks'] as const,
    storage: ['stats', 'storage'] as const,
    trends: (metric: string, period?: string) => ['stats', 'trends', metric, period] as const,
  },
  // Analytics
  analytics: {
    dashboard: ['analytics', 'dashboard'] as const,
    timeSeries: (metric: string, range?: string) => ['analytics', 'timeseries', metric, range] as const,
    projectPerformance: (projectId: string) => ['analytics', 'project-performance', projectId] as const,
  },
  // Admin
  admin: {
    users: (params?: Record<string, any>) => ['admin', 'users', params] as const,
    user: (id: string) => ['admin', 'user', id] as const,
    config: ['admin', 'config'] as const,
    health: ['admin', 'health'] as const,
    metrics: ['admin', 'metrics'] as const,
    alerts: (resolved?: boolean) => ['admin', 'alerts', resolved] as const,
    auditLogs: (params?: Record<string, any>) => ['admin', 'audit-logs', params] as const,
    backups: ['admin', 'backups'] as const,
    jobs: (status?: string) => ['admin', 'jobs', status] as const,
    resources: ['admin', 'resources'] as const,
  },
  // AI
  ai: {
    status: ['ai', 'status'] as const,
    conversations: ['ai', 'conversations'] as const,
    conversation: (id: string) => ['ai', 'conversation', id] as const,
    insights: (projectId: string) => ['ai', 'insights', projectId] as const,
  },
}
