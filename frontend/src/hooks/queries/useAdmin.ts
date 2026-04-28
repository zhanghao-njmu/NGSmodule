/**
 * TanStack Query hooks for the Admin domain.
 *
 * Covers user management, system config, audit logs, alerts, backups,
 * and the enhanced admin endpoints.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import adminService from '@/services/admin.service'
import enhancedAdminService from '@/services/admin.enhanced.service'

// ----- users --------------------------------------------------------------

export function useAdminUsers(params?: any) {
  return useQuery({
    queryKey: queryKeys.admin.users(params),
    queryFn: () => adminService.getUsers(params),
  })
}

export function useAdminUser(id: string | undefined) {
  return useQuery({
    queryKey: queryKeys.admin.user(id ?? ''),
    queryFn: () => adminService.getUser(id as string),
    enabled: !!id,
  })
}

export function useUpdateAdminUser() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, data }: { id: string; data: any }) => adminService.updateUser(id, data),
    onSuccess: (_, { id }) => {
      queryClient.invalidateQueries({ queryKey: queryKeys.admin.user(id) })
      queryClient.invalidateQueries({ queryKey: ['admin', 'users'] })
    },
  })
}

export function useChangeUserRole() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, role }: { id: string; role: 'user' | 'admin' }) => adminService.changeUserRole(id, role),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'users'] })
    },
  })
}

export function useToggleUserStatus() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, isActive, reason }: { id: string; isActive: boolean; reason?: string }) =>
      adminService.toggleUserStatus(id, isActive, reason),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'users'] })
    },
  })
}

export function useResetUserPassword() {
  return useMutation({
    mutationFn: ({ id, newPassword, notifyUser }: { id: string; newPassword: string; notifyUser?: boolean }) =>
      adminService.resetUserPassword(id, newPassword, notifyUser),
  })
}

export function useDeleteAdminUser() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: ({ id, options }: { id: string; options?: { transferDataTo?: string; reason?: string } }) =>
      adminService.deleteUser(id, options),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'users'] })
    },
  })
}

// ----- system config ------------------------------------------------------

export function useSystemConfig() {
  return useQuery({
    queryKey: queryKeys.admin.config,
    queryFn: () => adminService.getConfig(),
    staleTime: 60_000,
  })
}

export function useUpdateSystemConfig() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => adminService.updateConfig(data),
    onSuccess: (data) => {
      queryClient.setQueryData(queryKeys.admin.config, data)
    },
  })
}

export function useResetSystemConfig() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (categories?: string[]) => adminService.resetConfig(categories),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.admin.config })
    },
  })
}

// ----- system health & stats ----------------------------------------------

export function useSystemHealth() {
  return useQuery({
    queryKey: queryKeys.admin.health,
    queryFn: () => adminService.getSystemHealth(),
    staleTime: 30_000,
    refetchInterval: 60_000,
  })
}

export function useSystemStats() {
  return useQuery({
    queryKey: ['admin', 'system-stats'],
    queryFn: () => adminService.getSystemStats(),
    staleTime: 60_000,
  })
}

export function useAdminLogs(params?: any) {
  return useQuery({
    queryKey: ['admin', 'logs', params],
    queryFn: () => adminService.getLogs(params),
    staleTime: 30_000,
  })
}

export function useCleanupSystem() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (options: any) => adminService.cleanupSystem(options),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin'] })
    },
  })
}

// ----- enhanced: alerts / audit / backups / jobs --------------------------

export function useSystemMetrics() {
  return useQuery({
    queryKey: queryKeys.admin.metrics,
    queryFn: () => enhancedAdminService.getSystemMetrics(),
    staleTime: 15_000,
    refetchInterval: 30_000,
  })
}

export function useAlerts(resolved = false) {
  return useQuery({
    queryKey: queryKeys.admin.alerts(resolved),
    queryFn: () => enhancedAdminService.getAlerts(resolved),
    staleTime: 15_000,
    refetchInterval: 60_000,
  })
}

export function useResolveAlert() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (alertId: string) => enhancedAdminService.resolveAlert(alertId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'alerts'] })
    },
  })
}

export function useAuditLogs(filters?: any) {
  return useQuery({
    queryKey: queryKeys.admin.auditLogs(filters),
    queryFn: () => enhancedAdminService.getAuditLogs(filters),
  })
}

export function useExportAuditLogs() {
  return useMutation({
    mutationFn: (filters?: any) => enhancedAdminService.exportAuditLogs(filters),
  })
}

export function useResourceUsage() {
  return useQuery({
    queryKey: queryKeys.admin.resources,
    queryFn: () => enhancedAdminService.getResourceUsage(),
    staleTime: 30_000,
    refetchInterval: 60_000,
  })
}

export function useBackups() {
  return useQuery({
    queryKey: queryKeys.admin.backups,
    queryFn: () => enhancedAdminService.listBackups(),
  })
}

export function useCreateBackup() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (type: 'full' | 'incremental') => enhancedAdminService.createBackup(type),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.admin.backups })
    },
  })
}

export function useAdminJobs(status?: string) {
  return useQuery({
    queryKey: queryKeys.admin.jobs(status),
    queryFn: () => enhancedAdminService.listJobs(status),
    staleTime: 15_000,
    refetchInterval: 30_000,
  })
}

export function useCancelAdminJob() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (jobId: string) => enhancedAdminService.cancelJob(jobId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'jobs'] })
    },
  })
}

export function useRetryAdminJob() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (jobId: string) => enhancedAdminService.retryJob(jobId),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['admin', 'jobs'] })
    },
  })
}
