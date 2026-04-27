/**
 * TanStack Query hooks for the Notifications domain.
 *
 * Pairs with `useNotificationStream` (in `useRealtime.ts`) which invalidates
 * these caches whenever the backend pushes a new notification via WebSocket.
 */
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'

import { queryKeys } from '@/lib/queryClient'
import notificationService from '@/services/notification.service'
import type { NotificationListParams } from '@/types/notification'

export function useNotifications(params?: NotificationListParams) {
  return useQuery({
    queryKey: queryKeys.notifications.list(params as any),
    queryFn: () => notificationService.getNotifications(params),
    staleTime: 15_000,
  })
}

export function useUnreadNotificationCount() {
  return useQuery({
    queryKey: queryKeys.notifications.unreadCount,
    queryFn: () => notificationService.getUnreadCount(),
    staleTime: 10_000,
    // Light polling as a safety net even if WebSocket is down
    refetchInterval: 60_000,
  })
}

export function useNotificationSettings() {
  return useQuery({
    queryKey: queryKeys.notifications.settings,
    queryFn: () => notificationService.getSettings(),
    staleTime: 5 * 60_000,
  })
}

// ----- mutations -----

export function useMarkNotificationRead() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => notificationService.markAsRead(id),
    // Optimistic update — flip the read flag immediately so the UI
    // doesn't lag for a second.
    onMutate: async (id) => {
      await queryClient.cancelQueries({ queryKey: queryKeys.notifications.all })
      const previous = queryClient.getQueriesData({ queryKey: queryKeys.notifications.all })

      queryClient.setQueriesData<any>({ queryKey: queryKeys.notifications.all }, (old: any) => {
        if (!old?.items) {
          return old
        }
        return {
          ...old,
          items: old.items.map((n: any) => (n.id === id ? { ...n, read: true } : n)),
          unread_count: Math.max(0, (old.unread_count ?? 0) - 1),
        }
      })

      return { previous }
    },
    onError: (_err, _id, context: any) => {
      // Rollback
      if (context?.previous) {
        for (const [key, value] of context.previous) {
          queryClient.setQueryData(key, value)
        }
      }
    },
    onSettled: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.unreadCount })
    },
  })
}

export function useMarkAllNotificationsRead() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: () => notificationService.markAllAsRead(),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.unreadCount })
    },
  })
}

export function useDeleteNotification() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (id: string) => notificationService.deleteNotification(id),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.all })
      queryClient.invalidateQueries({ queryKey: queryKeys.notifications.unreadCount })
    },
  })
}

export function useUpdateNotificationSettings() {
  const queryClient = useQueryClient()
  return useMutation({
    mutationFn: (data: any) => notificationService.updateSettings(data),
    onSuccess: (data) => {
      queryClient.setQueryData(queryKeys.notifications.settings, data)
    },
  })
}
