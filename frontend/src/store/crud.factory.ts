/**
 * CRUD Store Factory
 * Generic factory for creating Zustand stores with CRUD operations
 * Eliminates code duplication across state management layers
 */

import { create } from 'zustand'
import { message } from 'antd'
import { CrudService } from '@/services/crud.factory'
import { ListParams } from '@/types/common'

/**
 * Generic CRUD store state
 * Common state properties for all CRUD stores
 */
export interface CrudStoreState<T> {
  /** Array of items */
  items: T[]

  /** Currently selected/viewed item */
  current: T | null

  /** Loading state */
  loading: boolean

  /** Error message */
  error: string | null

  /** Total number of items (for pagination) */
  total: number

  /** Pagination metadata */
  pagination: {
    page: number
    pageSize: number
  }
}

/**
 * Generic CRUD store actions
 * Common actions for all CRUD stores
 */
export interface CrudStoreActions<T, CreateT = Partial<T>, UpdateT = Partial<T>> {
  /**
   * Fetch all items with optional filters/pagination
   */
  fetchItems: (params?: ListParams) => Promise<void>

  /**
   * Fetch a single item by ID
   */
  fetchById: (id: string) => Promise<void>

  /**
   * Create a new item
   */
  createItem: (data: CreateT) => Promise<void>

  /**
   * Update an existing item
   */
  updateItem: (id: string, data: UpdateT) => Promise<void>

  /**
   * Delete an item by ID
   */
  deleteItem: (id: string) => Promise<void>

  /**
   * Set the current item
   */
  setCurrent: (item: T | null) => void

  /**
   * Clear error state
   */
  clearError: () => void

  /**
   * Reset store to initial state
   */
  reset: () => void

  /**
   * Set loading state
   */
  setLoading: (loading: boolean) => void

  /**
   * Set pagination
   */
  setPagination: (page: number, pageSize: number) => void
}

/**
 * Complete CRUD store type
 */
export type CrudStore<T, CreateT = Partial<T>, UpdateT = Partial<T>> = CrudStoreState<T> &
  CrudStoreActions<T, CreateT, UpdateT>

/**
 * Configuration options for CRUD store factory
 */
export interface CrudStoreConfig {
  /**
   * Human-readable entity name for success/error messages
   * @example 'project', 'sample', 'task'
   */
  entityName: string

  /**
   * Plural form of entity name for messages
   * @example 'projects', 'samples', 'tasks'
   * If not provided, will append 's' to entityName
   */
  entityNamePlural?: string

  /**
   * Whether to show success messages for CRUD operations
   * @default true
   */
  showSuccessMessages?: boolean

  /**
   * Whether to show error messages for CRUD operations
   * @default true
   */
  showErrorMessages?: boolean

  /**
   * Custom success messages
   */
  messages?: {
    created?: string
    updated?: string
    deleted?: string
    fetched?: string
  }

  /**
   * Initial pagination settings
   */
  initialPagination?: {
    page: number
    pageSize: number
  }
}

/**
 * Initial state factory
 */
const createInitialState = <T>(config?: CrudStoreConfig): CrudStoreState<T> => ({
  items: [],
  current: null,
  loading: false,
  error: null,
  total: 0,
  pagination: {
    page: config?.initialPagination?.page || 1,
    pageSize: config?.initialPagination?.pageSize || 10,
  },
})

/**
 * Create a CRUD store for a given entity type
 *
 * @param service - CRUD service instance
 * @param config - Store configuration
 * @returns Zustand store hook
 *
 * @example
 * // Create service
 * const projectService = createCrudService<Project>({ endpoint: 'projects' })
 *
 * // Create store
 * export const useProjectStore = createCrudStore(
 *   projectService,
 *   { entityName: 'project' }
 * )
 *
 * // Use in component
 * const { items, loading, fetchItems, createItem } = useProjectStore()
 */
export function createCrudStore<T extends { id: string }, CreateT = Partial<T>, UpdateT = Partial<T>>(
  service: CrudService<T, CreateT, UpdateT>,
  config: CrudStoreConfig
) {
  const {
    entityName,
    entityNamePlural = `${entityName}s`,
    showSuccessMessages = true,
    showErrorMessages = true,
    messages = {},
  } = config

  const initialState = createInitialState<T>(config)

  return create<CrudStore<T, CreateT, UpdateT>>((set, get) => ({
    // ========== State ==========
    ...initialState,

    // ========== Actions ==========

    /**
     * Fetch all items with optional filters and pagination
     */
    fetchItems: async (params?: ListParams) => {
      set({ loading: true, error: null })

      try {
        const response = await service.getAll(params)

        set({
          items: response.items,
          total: response.total,
          loading: false,
          error: null,
        })

        if (showSuccessMessages && messages.fetched) {
          message.success(messages.fetched)
        }
      } catch (error: any) {
        const errorMessage = error.message || `Failed to fetch ${entityNamePlural}`

        set({
          error: errorMessage,
          loading: false,
        })

        if (showErrorMessages) {
          message.error(errorMessage)
        }
      }
    },

    /**
     * Fetch a single item by ID
     */
    fetchById: async (id: string) => {
      set({ loading: true, error: null })

      try {
        const item = await service.getById(id)

        set({
          current: item,
          loading: false,
          error: null,
        })
      } catch (error: any) {
        const errorMessage = error.message || `Failed to fetch ${entityName}`

        set({
          error: errorMessage,
          loading: false,
        })

        if (showErrorMessages) {
          message.error(errorMessage)
        }
      }
    },

    /**
     * Create a new item
     */
    createItem: async (data: CreateT) => {
      set({ loading: true, error: null })

      try {
        const newItem = await service.create(data)

        set((state) => ({
          items: [newItem, ...state.items],
          total: state.total + 1,
          loading: false,
          error: null,
        }))

        const successMessage =
          messages.created || `${entityName.charAt(0).toUpperCase() + entityName.slice(1)} created successfully`

        if (showSuccessMessages) {
          message.success(successMessage)
        }
      } catch (error: any) {
        const errorMessage = error.message || `Failed to create ${entityName}`

        set({
          error: errorMessage,
          loading: false,
        })

        if (showErrorMessages) {
          message.error(errorMessage)
        }

        throw error // Re-throw to allow caller to handle
      }
    },

    /**
     * Update an existing item
     */
    updateItem: async (id: string, data: UpdateT) => {
      set({ loading: true, error: null })

      try {
        const updatedItem = await service.update(id, data)

        set((state) => ({
          items: state.items.map((item) => (item.id === id ? updatedItem : item)),
          current: state.current?.id === id ? updatedItem : state.current,
          loading: false,
          error: null,
        }))

        const successMessage =
          messages.updated || `${entityName.charAt(0).toUpperCase() + entityName.slice(1)} updated successfully`

        if (showSuccessMessages) {
          message.success(successMessage)
        }
      } catch (error: any) {
        const errorMessage = error.message || `Failed to update ${entityName}`

        set({
          error: errorMessage,
          loading: false,
        })

        if (showErrorMessages) {
          message.error(errorMessage)
        }

        throw error // Re-throw to allow caller to handle
      }
    },

    /**
     * Delete an item by ID
     */
    deleteItem: async (id: string) => {
      set({ loading: true, error: null })

      try {
        await service.delete(id)

        set((state) => ({
          items: state.items.filter((item) => item.id !== id),
          total: state.total - 1,
          current: state.current?.id === id ? null : state.current,
          loading: false,
          error: null,
        }))

        const successMessage =
          messages.deleted || `${entityName.charAt(0).toUpperCase() + entityName.slice(1)} deleted successfully`

        if (showSuccessMessages) {
          message.success(successMessage)
        }
      } catch (error: any) {
        const errorMessage = error.message || `Failed to delete ${entityName}`

        set({
          error: errorMessage,
          loading: false,
        })

        if (showErrorMessages) {
          message.error(errorMessage)
        }

        throw error // Re-throw to allow caller to handle
      }
    },

    /**
     * Set the current item
     */
    setCurrent: (item: T | null) => {
      set({ current: item })
    },

    /**
     * Clear error state
     */
    clearError: () => {
      set({ error: null })
    },

    /**
     * Reset store to initial state
     */
    reset: () => {
      set(initialState)
    },

    /**
     * Set loading state
     */
    setLoading: (loading: boolean) => {
      set({ loading })
    },

    /**
     * Set pagination
     */
    setPagination: (page: number, pageSize: number) => {
      set((state) => ({
        pagination: { page, pageSize },
      }))
    },
  }))
}

/**
 * Utility: Extend a CRUD store with custom state and actions
 *
 * @param baseStore - Base CRUD store
 * @param customState - Additional state properties
 * @param customActions - Additional actions
 * @returns Extended store
 *
 * @example
 * const useProjectStore = extendCrudStore(
 *   baseCrudStore,
 *   { filters: { status: 'active' } },
 *   (set, get) => ({
 *     setFilter: (status: string) => set({ filters: { status } }),
 *     async getStats(projectId: string) {
 *       const stats = await projectService.getStats(projectId)
 *       set({ stats })
 *     }
 *   })
 * )
 */
export function extendCrudStore<
  BaseState,
  BaseActions,
  CustomState extends Record<string, any>,
  CustomActions extends Record<string, any>
>(
  baseStore: () => BaseState & BaseActions,
  customState: CustomState,
  customActions: (
    set: (partial: Partial<BaseState & CustomState>) => void,
    get: () => BaseState & BaseActions & CustomState & CustomActions
  ) => CustomActions
) {
  return create<BaseState & BaseActions & CustomState & CustomActions>((set, get) => ({
    ...(baseStore() as any),
    ...customState,
    ...customActions(set as any, get),
  }))
}

/**
 * Helper: Create a simple store (no service integration)
 * Useful for UI-only stores that don't need API integration
 *
 * @param initialState - Initial state
 * @param actions - Store actions
 * @returns Zustand store hook
 *
 * @example
 * export const useUIStore = createSimpleStore(
 *   { sidebarOpen: true, theme: 'light' },
 *   (set) => ({
 *     toggleSidebar: () => set((state) => ({ sidebarOpen: !state.sidebarOpen })),
 *     setTheme: (theme: string) => set({ theme })
 *   })
 * )
 */
export function createSimpleStore<State extends Record<string, any>, Actions extends Record<string, any>>(
  initialState: State,
  actions: (set: (partial: Partial<State>) => void, get: () => State & Actions) => Actions
) {
  return create<State & Actions>((set, get) => ({
    ...initialState,
    ...actions(set, get as () => State & Actions),
  }))
}

export default createCrudStore
