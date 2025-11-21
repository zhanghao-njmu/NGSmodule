/**
 * WebSocket Service - Real-time task updates
 */
import type { WebSocketMessage } from '../types/task'

type MessageHandler = (message: WebSocketMessage) => void

class WebSocketService {
  private ws: WebSocket | null = null
  private messageHandlers: Set<MessageHandler> = new Set()
  private reconnectAttempts = 0
  private maxReconnectAttempts = 5
  private reconnectDelay = 2000
  private pingInterval: number | null = null
  private url: string = ''

  /**
   * Connect to WebSocket server
   */
  connect(token: string): void {
    const wsUrl =
      import.meta.env.VITE_WS_URL ||
      `ws://${window.location.hostname}:8000/api/v1/ws?token=${token}`

    this.url = wsUrl
    this.ws = new WebSocket(wsUrl)

    this.ws.onopen = () => {
      console.log('✅ WebSocket connected')
      this.reconnectAttempts = 0

      // Start ping interval to keep connection alive
      this.startPing()
    }

    this.ws.onmessage = (event) => {
      try {
        const message: WebSocketMessage = JSON.parse(event.data)
        this.notifyHandlers(message)
      } catch (error) {
        console.error('Error parsing WebSocket message:', error)
      }
    }

    this.ws.onerror = (error) => {
      console.error('❌ WebSocket error:', error)
    }

    this.ws.onclose = () => {
      console.log('👋 WebSocket disconnected')
      this.stopPing()

      // Attempt to reconnect
      if (this.reconnectAttempts < this.maxReconnectAttempts) {
        this.reconnectAttempts++
        console.log(`🔄 Reconnecting... (attempt ${this.reconnectAttempts})`)
        setTimeout(() => {
          if (token) {
            this.connect(token)
          }
        }, this.reconnectDelay * this.reconnectAttempts)
      }
    }
  }

  /**
   * Disconnect from WebSocket server
   */
  disconnect(): void {
    this.stopPing()
    if (this.ws) {
      this.ws.close()
      this.ws = null
    }
  }

  /**
   * Subscribe to task updates
   */
  subscribeToTask(taskId: string): void {
    this.send({
      type: 'subscribe',
      task_id: taskId,
    })
  }

  /**
   * Unsubscribe from task updates
   */
  unsubscribeFromTask(taskId: string): void {
    this.send({
      type: 'unsubscribe',
      task_id: taskId,
    })
  }

  /**
   * Send ping to keep connection alive
   */
  private ping(): void {
    this.send({ type: 'ping' })
  }

  /**
   * Start ping interval
   */
  private startPing(): void {
    this.pingInterval = window.setInterval(() => {
      this.ping()
    }, 30000) // Ping every 30 seconds
  }

  /**
   * Stop ping interval
   */
  private stopPing(): void {
    if (this.pingInterval !== null) {
      clearInterval(this.pingInterval)
      this.pingInterval = null
    }
  }

  /**
   * Send message to WebSocket server
   */
  private send(message: Partial<WebSocketMessage>): void {
    if (this.ws && this.ws.readyState === WebSocket.OPEN) {
      this.ws.send(JSON.stringify(message))
    } else {
      console.warn('WebSocket is not connected')
    }
  }

  /**
   * Add message handler
   */
  addMessageHandler(handler: MessageHandler): void {
    this.messageHandlers.add(handler)
  }

  /**
   * Remove message handler
   */
  removeMessageHandler(handler: MessageHandler): void {
    this.messageHandlers.delete(handler)
  }

  /**
   * Notify all message handlers
   */
  private notifyHandlers(message: WebSocketMessage): void {
    this.messageHandlers.forEach((handler) => {
      try {
        handler(message)
      } catch (error) {
        console.error('Error in message handler:', error)
      }
    })
  }

  /**
   * Check if WebSocket is connected
   */
  isConnected(): boolean {
    return this.ws !== null && this.ws.readyState === WebSocket.OPEN
  }
}

export const websocketService = new WebSocketService()
export default websocketService
