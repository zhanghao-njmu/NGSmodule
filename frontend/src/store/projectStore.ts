/**
 * Project Store - Global state management for projects
 */
import { create } from 'zustand'
import { projectService } from '../services/project.service'
import type { Project, ProjectStats, ProjectCreate, ProjectUpdate } from '../types/project'
import { message } from 'antd'

interface ProjectStore {
  // State
  projects: Project[]
  currentProject: Project | null
  stats: ProjectStats | null
  loading: boolean
  error: string | null

  // Actions
  fetchProjects: (params?: { status?: string }) => Promise<void>
  fetchProjectById: (id: string) => Promise<void>
  fetchStats: () => Promise<void>
  createProject: (data: ProjectCreate) => Promise<Project | null>
  updateProject: (id: string, data: ProjectUpdate) => Promise<void>
  deleteProject: (id: string) => Promise<void>
  archiveProject: (id: string) => Promise<void>
  restoreProject: (id: string) => Promise<void>
  setCurrentProject: (project: Project | null) => void
  clearError: () => void
}

export const useProjectStore = create<ProjectStore>((set, get) => ({
  // Initial state
  projects: [],
  currentProject: null,
  stats: null,
  loading: false,
  error: null,

  // Fetch all projects
  fetchProjects: async (params) => {
    set({ loading: true, error: null })
    try {
      const response = await projectService.getProjects(params)
      set({ projects: response.items, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch projects: ${error.message}`)
    }
  },

  // Fetch project by ID
  fetchProjectById: async (id) => {
    set({ loading: true, error: null })
    try {
      const project = await projectService.getProject(id)
      set({ currentProject: project, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch project: ${error.message}`)
    }
  },

  // Fetch project statistics
  fetchStats: async () => {
    set({ loading: true, error: null })
    try {
      const stats = await projectService.getStats()
      set({ stats, loading: false })
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to fetch statistics: ${error.message}`)
    }
  },

  // Create new project
  createProject: async (data) => {
    set({ loading: true, error: null })
    try {
      const project = await projectService.createProject(data)
      set((state) => ({
        projects: [project, ...state.projects],
        loading: false,
      }))
      message.success('Project created successfully')
      return project
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to create project: ${error.message}`)
      return null
    }
  },

  // Update project
  updateProject: async (id, data) => {
    set({ loading: true, error: null })
    try {
      const updatedProject = await projectService.updateProject(id, data)
      set((state) => ({
        projects: state.projects.map((p) => (p.id === id ? updatedProject : p)),
        currentProject:
          state.currentProject?.id === id ? updatedProject : state.currentProject,
        loading: false,
      }))
      message.success('Project updated successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to update project: ${error.message}`)
    }
  },

  // Delete project
  deleteProject: async (id) => {
    set({ loading: true, error: null })
    try {
      await projectService.deleteProject(id)
      set((state) => ({
        projects: state.projects.filter((p) => p.id !== id),
        currentProject: state.currentProject?.id === id ? null : state.currentProject,
        loading: false,
      }))
      message.success('Project deleted successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to delete project: ${error.message}`)
    }
  },

  // Archive project
  archiveProject: async (id) => {
    set({ loading: true, error: null })
    try {
      await projectService.archiveProject(id)
      // Refresh projects list
      await get().fetchProjects()
      message.success('Project archived successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to archive project: ${error.message}`)
    }
  },

  // Restore project
  restoreProject: async (id) => {
    set({ loading: true, error: null })
    try {
      await projectService.restoreProject(id)
      // Refresh projects list
      await get().fetchProjects()
      message.success('Project restored successfully')
    } catch (error: any) {
      set({ error: error.message, loading: false })
      message.error(`Failed to restore project: ${error.message}`)
    }
  },

  // Set current project
  setCurrentProject: (project) => {
    set({ currentProject: project })
  },

  // Clear error
  clearError: () => {
    set({ error: null })
  },
}))

export default useProjectStore
