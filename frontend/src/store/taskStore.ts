/**
 * Task Store - Global state management for tasks
 */
import { create } from 'zustand'
import { toast } from '@/utils/notification'
import { taskService } from '../services/task.service'
import { websocketService } from '../services/websocket.service'
import type { Task, TaskStats, TaskCreate, TaskExecuteRequest, WebSocketMessage } from '../types/task'
import { message } from 'antd'

interface TaskStore {
  // State
  tasks: Task[]
  currentTask: Task | null
  stats: TaskStats | null
  loading: boolean
  error: string | null
  taskLogs: { [taskId: string]: string }

  // Actions
  fetchTasks: (params?: { project_id?: string; status?: string }) => Promise<void>
  fetchTaskById: (id: string) => Promise<void>
  fetchStats: (params?: { project_id?: string }) => Promise<void>
  createTask: (data: TaskCreate) => Promise<Task | null>
  executeTask: (id: string, data: TaskExecuteRequest) => Promise<void>
  cancelTask: (id: string) => Promise<void>
  deleteTask: (id: string) => Promise<void>
  fetchTaskLogs: (id: string) => Promise<void>
  subscribeToTask: (id: string) => void
  unsubscribeFromTask: (id: string) => void
  handleWebSocketMessage: (message: WebSocketMessage) => void
  setCurrentTask: (task: Task | null) => void
  clearError: () => void
}

export const useTaskStore = create<TaskStore>((set, get) => ({
  // Initial state
  tasks: [],
  currentTask: null,
  stats: null,
  loading: false,
  error: null,
  taskLogs: {},

  // Fetch all tasks
  fetchTasks: async (params) => {
    set({ loading: true, error: null })
    try {
      const response = await taskService.getAll(params)
      set({ tasks: response.items, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch tasks: ${error.message}`)
    }
  },

  // Fetch task by ID
  fetchTaskById: async (id) => {
    set({ loading: true, error: null })
    try {
      const task = await taskService.getById(id)
      set({ currentTask: task, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch task: ${error.message}`)
    }
  },

  // Fetch task statistics
  fetchStats: async (params) => {
    set({ loading: true, error: null })
    try {
      const stats = await taskService.getStats(params)
      set({ stats, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch statistics: ${error.message}`)
    }
  },

  // Create new task
  createTask: async (data) => {
    set({ loading: true, error: null })
    try {
      const task = await taskService.create(data)
      set((state) => ({
        tasks: [task, ...state.tasks],
        loading: false,
      }))
      message.success('Task created successfully')
      return task
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to create task: ${error.message}`)
      return null
    }
  },

  // Execute task
  executeTask: async (id, data) => {
    set({ loading: true, error: null })
    try {
      await taskService.executeTask(id, data)
      // Subscribe to task updates
      get().subscribeToTask(id)
      // Refresh task details
      await get().fetchTaskById(id)
      message.success('Task execution started')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to execute task: ${error.message}`)
    }
  },

  // Cancel task
  cancelTask: async (id) => {
    set({ loading: true, error: null })
    try {
      await taskService.cancelTask(id)
      // Refresh task details
      await get().fetchTaskById(id)
      message.success('Task cancelled successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to cancel task: ${error.message}`)
    }
  },

  // Delete task
  deleteTask: async (id) => {
    set({ loading: true, error: null })
    try {
      await taskService.delete(id)
      set((state) => ({
        tasks: state.tasks.filter((t) => t.id !== id),
        currentTask: state.currentTask?.id === id ? null : state.currentTask,
        loading: false,
      }))
      message.success('Task deleted successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to delete task: ${error.message}`)
    }
  },

  // Fetch task logs
  fetchTaskLogs: async (id) => {
    set({ loading: true, error: null })
    try {
      const response = await taskService.getTaskLogs(id)
      set((state) => ({
        taskLogs: {
          ...state.taskLogs,
          [id]: response.log_content,
        },
        loading: false,
      }))
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch task logs: ${error.message}`)
    }
  },

  // Subscribe to task updates via WebSocket
  subscribeToTask: (id) => {
    if (websocketService.isConnected()) {
      websocketService.subscribeToTask(id)
    }
  },

  // Unsubscribe from task updates
  unsubscribeFromTask: (id) => {
    if (websocketService.isConnected()) {
      websocketService.unsubscribeFromTask(id)
    }
  },

  // Handle WebSocket messages
  handleWebSocketMessage: (msg) => {
    if (msg.type === 'task_update' && msg.task_id) {
      // Update task in list
      set((state) => ({
        tasks: state.tasks.map((task) =>
          task.id === msg.task_id
            ? {
                ...task,
                status: msg.status || task.status,
                progress: msg.progress ?? task.progress,
              }
            : task,
        ),
        currentTask:
          state.currentTask?.id === msg.task_id
            ? {
                ...state.currentTask,
                status: msg.status || state.currentTask.status,
                progress: msg.progress ?? state.currentTask.progress,
              }
            : state.currentTask,
      }))

      // Show notification for completed/failed tasks
      if (msg.status === 'completed') {
        toast.success(`Task completed: ${msg.message}`)
      } else if (msg.status === 'failed') {
        toast.error(`Task failed: ${msg.message}`)
      }
    }
  },

  // Set current task
  setCurrentTask: (task) => {
    set({ currentTask: task })
  },

  // Clear error
  clearError: () => {
    set({ error: null })
  },
}))

export default useTaskStore
