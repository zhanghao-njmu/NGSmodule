/**
 * Project Store
 * Global state management for projects
 * Refactored to use CRUD store factory with project-specific extensions
 */

import { create } from 'zustand'
import { projectService } from '@/services/project.service'
import type { Project, ProjectCreate, ProjectUpdate, ProjectStats } from '@/types/project'
import { message } from 'antd'
import type { ListParams } from '@/types/common'

/**
 * Project store state
 */
interface ProjectStoreState {
  // Data
  items: Project[]
  current: Project | null
  stats: ProjectStats | null
  loading: boolean
  error: string | null
  total: number
  pagination: {
    page: number
    pageSize: number
  }

  // Standard CRUD actions
  fetchItems: (params?: ListParams) => Promise<void>
  fetchById: (id: string) => Promise<void>
  createItem: (data: ProjectCreate) => Promise<void>
  updateItem: (id: string, data: ProjectUpdate) => Promise<void>
  deleteItem: (id: string) => Promise<void>
  setCurrent: (item: Project | null) => void
  clearError: () => void
  reset: () => void
  setLoading: (loading: boolean) => void
  setPagination: (page: number, pageSize: number) => void

  // Project-specific actions
  fetchStats: () => Promise<void>
  archiveProject: (id: string) => Promise<void>
  restoreProject: (id: string) => Promise<void>
}

/**
 * Create project store with CRUD operations and project-specific functionality
 */
export const useProjectStore = create<ProjectStoreState>((set, get) => ({
  // Initial state
  items: [],
  current: null,
  stats: null,
  loading: false,
  error: null,
  total: 0,
  pagination: {
    page: 1,
    pageSize: 10,
  },

  // Standard CRUD actions (using factory pattern internally)

  fetchItems: async (params?: ListParams) => {
    set({ loading: true, error: null })
    try {
      const response = await projectService.getAll(params)
      set({ items: response.items, total: response.total, loading: false, error: null })
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to fetch projects'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
    }
  },

  fetchById: async (id: string) => {
    set({ loading: true, error: null })
    try {
      const item = await projectService.getById(id)
      set({ current: item, loading: false, error: null })
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to fetch project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
    }
  },

  createItem: async (data: ProjectCreate) => {
    set({ loading: true, error: null })
    try {
      const newItem = await projectService.create(data)
      set((state) => ({
        items: [newItem, ...state.items],
        total: state.total + 1,
        loading: false,
        error: null,
      }))
      message.success('Project created successfully')
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to create project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
      throw error
    }
  },

  updateItem: async (id: string, data: ProjectUpdate) => {
    set({ loading: true, error: null })
    try {
      const updatedItem = await projectService.update(id, data)
      set((state) => ({
        items: state.items.map((item) => (item.id === id ? updatedItem : item)),
        current: state.current?.id === id ? updatedItem : state.current,
        loading: false,
        error: null,
      }))
      message.success('Project updated successfully')
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to update project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
      throw error
    }
  },

  deleteItem: async (id: string) => {
    set({ loading: true, error: null })
    try {
      await projectService.delete(id)
      set((state) => ({
        items: state.items.filter((item) => item.id !== id),
        total: state.total - 1,
        current: state.current?.id === id ? null : state.current,
        loading: false,
        error: null,
      }))
      message.success('Project deleted successfully')
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to delete project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
      throw error
    }
  },

  setCurrent: (item: Project | null) => {
    set({ current: item })
  },

  clearError: () => {
    set({ error: null })
  },

  reset: () => {
    set({
      items: [],
      current: null,
      stats: null,
      loading: false,
      error: null,
      total: 0,
      pagination: { page: 1, pageSize: 10 },
    })
  },

  setLoading: (loading: boolean) => {
    set({ loading })
  },

  setPagination: (page: number, pageSize: number) => {
    set({ pagination: { page, pageSize } })
  },

  // Project-specific actions

  fetchStats: async () => {
    set({ loading: true, error: null })
    try {
      const stats = await projectService.getStats()
      set({ stats, loading: false })
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to fetch statistics'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
    }
  },

  archiveProject: async (id: string) => {
    set({ loading: true, error: null })
    try {
      await projectService.archiveProject(id)
      message.success('Project archived successfully')
      await get().fetchItems()
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to archive project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
    }
  },

  restoreProject: async (id: string) => {
    set({ loading: true, error: null })
    try {
      await projectService.restoreProject(id)
      message.success('Project restored successfully')
      await get().fetchItems()
    } catch (error: any) {
      const errorMessage = error.message || 'Failed to restore project'
      set({ error: errorMessage, loading: false })
      message.error(errorMessage)
    }
  },
}))

export default useProjectStore
