/**
 * Realtime WebSocket abstraction.
 *
 * Backend (P2-2) publishes events via Redis pub/sub:
 *   - realtime:user:{user_id}  -> personal events (notifications, alerts)
 *   - realtime:task:{task_id}  -> task progress
 *
 * The FastAPI WebSocket endpoint at /api/v1/ws bridges these to the
 * connected client. This module manages a single shared connection per
 * browser tab with auto-reconnect and a typed event bus.
 */

export type RealtimeEvent =
  | { type: 'task_update'; task_id: string; status: string; progress: number; message?: string }
  | { type: 'notification'; data: any }
  | { type: 'alert'; data: any }
  | { type: 'subscribed'; task_id: string; status?: string; progress?: number }
  | { type: 'unsubscribed'; task_id: string }
  | { type: 'pong'; timestamp: string }
  | { type: 'error'; message: string }
  | { type: string; [key: string]: any }

type Listener = (event: RealtimeEvent) => void

class RealtimeClient {
  private socket: WebSocket | null = null
  private listeners = new Set<Listener>()
  private reconnectAttempts = 0
  private reconnectTimer: number | null = null
  private heartbeatTimer: number | null = null
  private intentionallyClosed = false
  private token: string | null = null
  private subscribedTasks = new Set<string>()

  /**
   * Connect (or reconnect) using a JWT token. Called by the auth store
   * after login and on token refresh.
   */
  connect(token: string): void {
    this.token = token
    this.intentionallyClosed = false
    this.openSocket()
  }

  disconnect(): void {
    this.intentionallyClosed = true
    if (this.reconnectTimer) {
      window.clearTimeout(this.reconnectTimer)
      this.reconnectTimer = null
    }
    if (this.heartbeatTimer) {
      window.clearInterval(this.heartbeatTimer)
      this.heartbeatTimer = null
    }
    if (this.socket) {
      this.socket.close(1000, 'Client disconnect')
      this.socket = null
    }
    this.subscribedTasks.clear()
  }

  isConnected(): boolean {
    return this.socket?.readyState === WebSocket.OPEN
  }

  /** Subscribe a callback to all incoming events. Returns unsubscribe fn. */
  on(listener: Listener): () => void {
    this.listeners.add(listener)
    return () => this.listeners.delete(listener)
  }

  /** Subscribe to task progress updates for a specific task. */
  subscribeToTask(taskId: string): void {
    this.subscribedTasks.add(taskId)
    this.send({ type: 'subscribe', task_id: taskId })
  }

  unsubscribeFromTask(taskId: string): void {
    this.subscribedTasks.delete(taskId)
    this.send({ type: 'unsubscribe', task_id: taskId })
  }

  // ----- internals -----

  private openSocket(): void {
    if (!this.token) {
      return
    }

    // Build ws:// or wss:// URL based on current page protocol
    const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
    const host = window.location.host
    const url = `${protocol}//${host}/api/v1/ws?token=${encodeURIComponent(this.token)}`

    try {
      const socket = new WebSocket(url)
      this.socket = socket

      socket.onopen = () => {
        this.reconnectAttempts = 0
        this.startHeartbeat()
        // Re-subscribe to any tasks that were active before the disconnect
        for (const taskId of this.subscribedTasks) {
          socket.send(JSON.stringify({ type: 'subscribe', task_id: taskId }))
        }
      }

      socket.onmessage = (msg) => {
        try {
          const event: RealtimeEvent = JSON.parse(msg.data)
          // Backend may wrap in {type, data} envelope; normalize
          this.dispatch(event)
        } catch {
          // ignore malformed payloads
        }
      }

      socket.onerror = () => {
        // Errors trigger onclose right after; let the close handler reconnect
      }

      socket.onclose = () => {
        this.stopHeartbeat()
        this.socket = null
        if (!this.intentionallyClosed) {
          this.scheduleReconnect()
        }
      }
    } catch {
      this.scheduleReconnect()
    }
  }

  private scheduleReconnect(): void {
    if (this.reconnectTimer) {
      return
    }
    // Exponential backoff capped at 30s
    const delay = Math.min(1000 * Math.pow(2, this.reconnectAttempts), 30_000)
    this.reconnectAttempts++
    this.reconnectTimer = window.setTimeout(() => {
      this.reconnectTimer = null
      this.openSocket()
    }, delay)
  }

  private startHeartbeat(): void {
    this.stopHeartbeat()
    this.heartbeatTimer = window.setInterval(() => {
      this.send({ type: 'ping' })
    }, 25_000)
  }

  private stopHeartbeat(): void {
    if (this.heartbeatTimer) {
      window.clearInterval(this.heartbeatTimer)
      this.heartbeatTimer = null
    }
  }

  private send(payload: object): boolean {
    if (this.socket?.readyState === WebSocket.OPEN) {
      this.socket.send(JSON.stringify(payload))
      return true
    }
    return false
  }

  private dispatch(event: RealtimeEvent): void {
    for (const listener of this.listeners) {
      try {
        listener(event)
      } catch (err) {
        // Don't let one listener crash the others
        console.error('[realtime] listener error:', err)
      }
    }
  }
}

export const realtimeClient = new RealtimeClient()
