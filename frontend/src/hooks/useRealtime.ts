/**
 * React hooks that bridge the WebSocket realtime layer with TanStack Query
 * cache invalidation. Use these in feature pages instead of touching
 * `realtimeClient` directly.
 */
import { useEffect, useRef } from 'react'
import { useQueryClient } from '@tanstack/react-query'
import { message, notification } from 'antd'

import { realtimeClient, type RealtimeEvent } from '@/lib/realtime'
import { queryKeys } from '@/lib/queryClient'

/**
 * Subscribe to all realtime events. Most callers should prefer the more
 * targeted hooks below, but this is useful for debugging and global
 * dashboards.
 */
export function useRealtimeEvents(handler: (event: RealtimeEvent) => void, deps: React.DependencyList = []): void {
  const handlerRef = useRef(handler)
  useEffect(() => {
    handlerRef.current = handler
  }, [handler])

  useEffect(() => {
    return realtimeClient.on((event) => handlerRef.current(event))
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, deps)
}

/**
 * Listen for new notifications and invalidate the notifications cache so
 * the UI re-renders. Also pops a transient toast for high-priority
 * notifications.
 */
export function useNotificationStream(options: { showToast?: boolean } = { showToast: true }): void {
  const queryClient = useQueryClient()

  useRealtimeEvents((event) => {
    if (event.type !== 'notification') {
      return
    }

    // Invalidate caches that should reflect the new notification
    queryClient.invalidateQueries({ queryKey: queryKeys.notifications.all })
    queryClient.invalidateQueries({ queryKey: queryKeys.notifications.unreadCount })

    if (options.showToast) {
      const data = (event as any).data || {}
      const priority = data.priority || 'normal'
      const title = data.title || 'Notification'
      const body = data.message || ''

      if (priority === 'urgent' || priority === 'high') {
        notification.warning({
          message: title,
          description: body,
          placement: 'topRight',
          duration: 6,
        })
      } else {
        message.info({ content: title })
      }
    }
  })
}

/**
 * Subscribe to progress updates for a specific task. Updates the cached
 * task detail and the tasks list so all open views stay in sync.
 *
 * Pass `null` to unsubscribe (e.g., when the task page unmounts or the
 * user navigates away).
 */
export function useTaskProgress(taskId: string | null | undefined): void {
  const queryClient = useQueryClient()

  useEffect(() => {
    if (!taskId) {
      return
    }

    realtimeClient.subscribeToTask(taskId)
    return () => {
      realtimeClient.unsubscribeFromTask(taskId)
    }
  }, [taskId])

  useRealtimeEvents(
    (event) => {
      if (!taskId) {
        return
      }
      if (event.type !== 'task_update' || event.task_id !== taskId) {
        return
      }

      // Optimistically update the cached task detail
      queryClient.setQueryData<any>(queryKeys.tasks.detail(taskId), (old: any) => {
        if (!old) {
          return old
        }
        return {
          ...old,
          status: event.status,
          progress: event.progress,
          ...(event.message ? { message: event.message } : {}),
        }
      })

      // Mark the task list stale so it refetches in the background
      queryClient.invalidateQueries({ queryKey: queryKeys.tasks.all })

      // If the task finished, refresh results and stats too
      if (event.status === 'completed' || event.status === 'failed') {
        queryClient.invalidateQueries({ queryKey: queryKeys.results.all })
        queryClient.invalidateQueries({ queryKey: queryKeys.stats.all })
      }
    },
    [taskId],
  )
}

/**
 * Listen for vendor data download progress and refresh the data-downloads
 * jobs query whenever the worker pushes an update.
 *
 * The worker emits 'download_update' events on the per-task channel
 * (`realtime:task:<job_id>`); we subscribe to each job id whose status
 * is still in flight, then invalidate the cache on every event.
 */
export function useDownloadJobsRealtime(jobIds: string[]): void {
  const queryClient = useQueryClient()

  // Subscribe to each in-flight job's task channel.
  useEffect(() => {
    jobIds.forEach((id) => realtimeClient.subscribeToTask(id))
    return () => {
      jobIds.forEach((id) => realtimeClient.unsubscribeFromTask(id))
    }
    // Stringify for stable dep — same set of ids ⇒ same effect.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [jobIds.join(',')])

  useRealtimeEvents(
    (event) => {
      if (event.type !== 'download_update') {
        return
      }
      // Refresh the user's job list; React Query handles the diff.
      queryClient.invalidateQueries({ queryKey: ['data-downloads', 'jobs'] })
    },
    [queryClient],
  )
}

/**
 * Connect/disconnect lifecycle hook. Mount this once in the authenticated
 * shell layout so the WebSocket follows the user's session.
 */
export function useRealtimeConnection(token: string | null): void {
  useEffect(() => {
    if (!token) {
      realtimeClient.disconnect()
      return
    }
    realtimeClient.connect(token)
    return () => {
      // Don't disconnect on every re-render; only on token change/logout.
      // The cleanup runs only when the token actually changes.
    }
  }, [token])
}
