# Phase 22: Code Cleanup & Redundancy Removal Plan

## Executive Summary

This document outlines the comprehensive code redundancy cleanup plan for the NGSmodule frontend codebase. The analysis identified **15 redundancy issues** across 5 categories, with **4 high-severity** issues requiring immediate attention.

**Estimated Impact:**
- Code reduction: ~2,000+ lines of duplicate code
- Improved maintainability: Single source of truth for common patterns
- Better type safety: Generic types reduce errors
- Faster development: Reusable utilities and factories

---

## Redundancy Analysis Results

### High Severity Issues (4)
1. ✅ **Duplicate Global CSS Files** - Two global.css files with overlapping styles
2. ✅ **Repeated List Page Patterns** - 5 list pages with identical logic (~150 lines each)
3. ✅ **Identical Service Structure** - 5 services with duplicate CRUD methods
4. ✅ **Identical Store Patterns** - Multiple stores with duplicate state management

### Medium Severity Issues (7)
5. ✅ **Duplicate Empty State Components** - EnhancedEmptyState vs EmptyState
6. ✅ **Repeated Table Columns** - Action columns duplicated across pages
7. ✅ **Duplicate CSS Variables** - Variables defined in multiple files
8. ✅ **Repeated ListResponse Types** - Each entity defines own ListResponse
9. ✅ **Duplicate File Size Formatter** - Function duplicated across components
10. ✅ **Repeated Form Validations** - Validation rules repeated in forms
11. ✅ **Inconsistent Service Pattern** - ResultService uses class, others use objects

### Low Severity Issues (4)
12. ✅ **Repeated Date Formatting** - dayjs patterns duplicated
13. ✅ **Repeated Stats Pattern** - Similar stats interfaces
14. ✅ **Repeated Create/Update Types** - Pattern not using TS utilities
15. ✅ **Utility Classes** - Could use design tokens more effectively

---

## Execution Plan

### Phase 1: Critical Infrastructure (Priority 1)

#### Task 1.1: Consolidate Global CSS Files
**Files to modify:**
- DELETE: `frontend/src/assets/styles/global.css`
- KEEP: `frontend/src/styles/global.css`
- UPDATE: All imports referencing old path

**Actions:**
1. Audit all imports of `assets/styles/global.css`
2. Merge unique styles into `styles/global.css`
3. Delete `assets/styles/global.css`
4. Update all import statements

**Impact:** Eliminates ~100 lines of duplicate CSS

---

#### Task 1.2: Consolidate CSS Variables
**Files to modify:**
- `frontend/src/styles/variables.css` (single source of truth)
- `frontend/src/assets/styles/global.css` (remove variables)

**Actions:**
1. Audit all CSS variable usage
2. Standardize variable names in `styles/variables.css`
3. Remove duplicate variables from other files
4. Document variable naming convention

**Impact:** Single source of truth for design tokens

---

#### Task 1.3: Create Generic Types
**New file:** `frontend/src/types/common.ts`

**Content:**
```typescript
// Generic paginated response
export interface PaginatedResponse<T> {
  total: number
  items: T[]
  page?: number
  page_size?: number
  skip?: number
  limit?: number
}

// Generic list params
export interface ListParams {
  page?: number
  page_size?: number
  skip?: number
  limit?: number
  search?: string
  sort_by?: string
  sort_order?: 'asc' | 'desc'
}

// Generic API response
export interface ApiResponse<T> {
  success: boolean
  data: T
  message?: string
}

// Generic error response
export interface ApiError {
  message: string
  code?: string
  details?: Record<string, any>
}
```

**Files to update:**
- `frontend/src/types/project.ts` - Replace ProjectListResponse
- `frontend/src/types/sample.ts` - Replace SampleListResponse
- `frontend/src/types/task.ts` - Replace TaskListResponse
- `frontend/src/types/file.ts` - Replace FileListResponse
- `frontend/src/types/result.ts` - Replace ResultListResponse

**Impact:** Reduces type definitions by ~50 lines, improves consistency

---

#### Task 1.4: Create Formatting Utilities
**New file:** `frontend/src/utils/format.ts`

**Content:**
```typescript
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(relativeTime)

/**
 * Format file size in bytes to human-readable string
 */
export const formatFileSize = (bytes: number): string => {
  if (bytes === 0) return '0 B'
  const k = 1024
  const sizes = ['B', 'KB', 'MB', 'GB', 'TB']
  const i = Math.floor(Math.log(bytes) / Math.log(k))
  return `${(bytes / Math.pow(k, i)).toFixed(2)} ${sizes[i]}`
}

/**
 * Format date to standard format (YYYY-MM-DD HH:mm)
 */
export const formatDateTime = (date: string | Date): string => {
  return dayjs(date).format('YYYY-MM-DD HH:mm')
}

/**
 * Format date to relative time (e.g., "2 hours ago")
 */
export const formatDateRelative = (date: string | Date): string => {
  return dayjs(date).fromNow()
}

/**
 * Format date to long format (YYYY-MM-DD HH:mm:ss)
 */
export const formatDateLong = (date: string | Date): string => {
  return dayjs(date).format('YYYY-MM-DD HH:mm:ss')
}

/**
 * Format number with specified decimal places
 */
export const formatNumber = (value: number, decimals: number = 2): string => {
  return value.toFixed(decimals)
}

/**
 * Format percentage
 */
export const formatPercent = (value: number, decimals: number = 1): string => {
  return `${(value * 100).toFixed(decimals)}%`
}

/**
 * Format duration in seconds to human-readable string
 */
export const formatDuration = (seconds: number): string => {
  if (seconds < 60) return `${seconds}s`
  const minutes = Math.floor(seconds / 60)
  const remainingSeconds = seconds % 60
  if (minutes < 60) return `${minutes}m ${remainingSeconds}s`
  const hours = Math.floor(minutes / 60)
  const remainingMinutes = minutes % 60
  return `${hours}h ${remainingMinutes}m`
}
```

**Files to update:**
- `frontend/src/pages/files/FileList.tsx` - Use formatFileSize
- `frontend/src/pages/admin/AdminDashboard.tsx` - Use formatFileSize
- All list pages - Use date formatting utilities

**Impact:** Eliminates duplicate formatting code, ensures consistency

---

### Phase 2: Service Layer Optimization (Priority 2)

#### Task 2.1: Create Generic CRUD Service Factory
**New file:** `frontend/src/services/crud.factory.ts`

**Content:**
```typescript
import apiClient from './api'
import { PaginatedResponse, ListParams } from '@/types/common'

export interface CrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>> {
  getAll: (params?: ListParams) => Promise<PaginatedResponse<T>>
  getById: (id: string) => Promise<T>
  create: (data: CreateT) => Promise<T>
  update: (id: string, data: UpdateT) => Promise<T>
  delete: (id: string) => Promise<void>
}

/**
 * Create a CRUD service for a given endpoint
 */
export function createCrudService<T, CreateT = Partial<T>, UpdateT = Partial<T>>(
  endpoint: string
): CrudService<T, CreateT, UpdateT> {
  return {
    getAll: async (params?: ListParams) => {
      return apiClient.get<PaginatedResponse<T>>(`/${endpoint}`, { params })
    },

    getById: async (id: string) => {
      return apiClient.get<T>(`/${endpoint}/${id}`)
    },

    create: async (data: CreateT) => {
      return apiClient.post<T>(`/${endpoint}`, data)
    },

    update: async (id: string, data: UpdateT) => {
      return apiClient.put<T>(`/${endpoint}/${id}`, data)
    },

    delete: async (id: string) => {
      return apiClient.delete(`/${endpoint}/${id}`)
    },
  }
}
```

**Files to refactor:**
- `frontend/src/services/project.service.ts` - Use factory + extend
- `frontend/src/services/sample.service.ts` - Use factory + extend
- `frontend/src/services/task.service.ts` - Use factory + extend
- `frontend/src/services/file.service.ts` - Use factory + extend
- `frontend/src/services/pipeline.service.ts` - Use factory + extend

**Impact:** Reduces service code from ~75 lines to ~20 lines each

---

#### Task 2.2: Standardize ResultService
**File:** `frontend/src/services/result.service.ts`

**Actions:**
1. Convert from class-based to object literal pattern
2. Use createCrudService factory as base
3. Extend with domain-specific methods

**Impact:** Consistency across all services

---

#### Task 2.3: Create Validation Rules Library
**New file:** `frontend/src/utils/validationRules.ts`

**Content:**
```typescript
import type { Rule } from 'antd/es/form'

export const validationRules = {
  /**
   * Required field rule
   */
  required: (fieldName: string): Rule => ({
    required: true,
    message: `Please input ${fieldName}!`,
  }),

  /**
   * Email validation rule
   */
  email: {
    type: 'email' as const,
    message: 'Please input a valid email address!',
  },

  /**
   * Minimum length rule
   */
  minLength: (min: number): Rule => ({
    min,
    message: `Must be at least ${min} characters!`,
  }),

  /**
   * Maximum length rule
   */
  maxLength: (max: number): Rule => ({
    max,
    message: `Must be at most ${max} characters!`,
  }),

  /**
   * Username validation rules
   */
  username: [
    { required: true, message: 'Please input your username!' },
    { min: 3, message: 'Username must be at least 3 characters!' },
    { max: 50, message: 'Username must be at most 50 characters!' },
    { pattern: /^[a-zA-Z0-9_-]+$/, message: 'Username can only contain letters, numbers, hyphens, and underscores!' },
  ],

  /**
   * Password validation rules
   */
  password: [
    { required: true, message: 'Please input your password!' },
    { min: 8, message: 'Password must be at least 8 characters!' },
    { pattern: /[A-Z]/, message: 'Password must contain at least one uppercase letter!' },
    { pattern: /[a-z]/, message: 'Password must contain at least one lowercase letter!' },
    { pattern: /[0-9]/, message: 'Password must contain at least one number!' },
  ],

  /**
   * Confirm password rule
   */
  confirmPassword: (getFieldValue: (field: string) => any): Rule => ({
    validator: (_, value) => {
      if (!value || getFieldValue('password') === value) {
        return Promise.resolve()
      }
      return Promise.reject(new Error('The two passwords do not match!'))
    },
  }),

  /**
   * URL validation rule
   */
  url: {
    type: 'url' as const,
    message: 'Please input a valid URL!',
  },

  /**
   * Number range rule
   */
  numberRange: (min: number, max: number): Rule => ({
    type: 'number' as const,
    min,
    max,
    message: `Must be between ${min} and ${max}!`,
  }),
}
```

**Files to update:**
- `frontend/src/pages/auth/Register.tsx` - Use validation rules
- `frontend/src/pages/auth/Login.tsx` - Use validation rules
- All form components - Use validation rules

**Impact:** Reduces form validation code duplication

---

### Phase 3: Component Optimization (Priority 3)

#### Task 3.1: Create List Page Hooks
**New file:** `frontend/src/hooks/useListPage.ts`

**Content:**
```typescript
import { useState, useEffect } from 'react'
import { Form, message } from 'antd'
import { PaginatedResponse } from '@/types/common'

export interface UseListPageOptions<T> {
  fetchData: (params?: any) => Promise<PaginatedResponse<T>>
  deleteItem?: (id: string) => Promise<void>
  onError?: (error: Error) => void
}

export function useListPage<T extends { id: string }>({
  fetchData,
  deleteItem,
  onError,
}: UseListPageOptions<T>) {
  const [items, setItems] = useState<T[]>([])
  const [loading, setLoading] = useState(false)
  const [total, setTotal] = useState(0)
  const [page, setPage] = useState(1)
  const [pageSize, setPageSize] = useState(10)
  const [searchText, setSearchText] = useState('')

  const loadData = async (params?: any) => {
    setLoading(true)
    try {
      const response = await fetchData({
        page,
        page_size: pageSize,
        search: searchText,
        ...params,
      })
      setItems(response.items)
      setTotal(response.total)
    } catch (error: any) {
      message.error(`Failed to load data: ${error.message}`)
      onError?.(error)
    } finally {
      setLoading(false)
    }
  }

  useEffect(() => {
    loadData()
  }, [page, pageSize, searchText])

  const handleDelete = async (id: string) => {
    if (!deleteItem) return
    try {
      await deleteItem(id)
      message.success('Deleted successfully')
      loadData()
    } catch (error: any) {
      message.error(`Failed to delete: ${error.message}`)
    }
  }

  const handleSearch = (value: string) => {
    setSearchText(value)
    setPage(1)
  }

  const handlePageChange = (newPage: number, newPageSize?: number) => {
    setPage(newPage)
    if (newPageSize) setPageSize(newPageSize)
  }

  const refresh = () => loadData()

  return {
    items,
    loading,
    total,
    page,
    pageSize,
    searchText,
    handleDelete,
    handleSearch,
    handlePageChange,
    refresh,
  }
}
```

**New file:** `frontend/src/hooks/useModal.ts`

**Content:**
```typescript
import { useState } from 'react'
import { Form } from 'antd'

export function useModal<T>(options?: {
  onSubmit?: (values: T, editing?: T | null) => Promise<void>
  onCancel?: () => void
}) {
  const [visible, setVisible] = useState(false)
  const [editing, setEditing] = useState<T | null>(null)
  const [form] = Form.useForm()

  const open = (item?: T) => {
    if (item) {
      setEditing(item)
      form.setFieldsValue(item)
    } else {
      setEditing(null)
      form.resetFields()
    }
    setVisible(true)
  }

  const close = () => {
    setVisible(false)
    setEditing(null)
    form.resetFields()
    options?.onCancel?.()
  }

  const submit = async (values: T) => {
    await options?.onSubmit?.(values, editing)
    close()
  }

  return {
    visible,
    editing,
    form,
    open,
    close,
    submit,
  }
}
```

**Impact:** Each list page reduces from ~450 lines to ~200 lines

---

#### Task 3.2: Create Table Column Factories
**New file:** `frontend/src/utils/tableColumns.tsx`

**Content:**
```typescript
import { Space, Button, Tooltip, Popconfirm, Tag } from 'antd'
import { EditOutlined, DeleteOutlined, EyeOutlined } from '@ant-design/icons'
import type { ColumnType } from 'antd/es/table'
import { formatDateTime, formatDateRelative, formatFileSize } from './format'

export const createActionColumn = <T extends { id: string }>(
  onEdit?: (record: T) => void,
  onDelete?: (record: T) => void,
  onView?: (record: T) => void,
  options?: { width?: number; fixed?: 'left' | 'right' }
): ColumnType<T> => ({
  title: 'Actions',
  key: 'actions',
  width: options?.width || 120,
  fixed: options?.fixed || 'right',
  render: (_, record) => (
    <Space>
      {onView && (
        <Tooltip title="View">
          <Button
            type="text"
            size="small"
            icon={<EyeOutlined />}
            onClick={() => onView(record)}
          />
        </Tooltip>
      )}
      {onEdit && (
        <Tooltip title="Edit">
          <Button
            type="text"
            size="small"
            icon={<EditOutlined />}
            onClick={() => onEdit(record)}
          />
        </Tooltip>
      )}
      {onDelete && (
        <Popconfirm
          title="Are you sure you want to delete this item?"
          onConfirm={() => onDelete(record)}
          okText="Yes"
          cancelText="No"
        >
          <Tooltip title="Delete">
            <Button type="text" size="small" danger icon={<DeleteOutlined />} />
          </Tooltip>
        </Popconfirm>
      )}
    </Space>
  ),
})

export const createDateColumn = <T,>(
  title: string,
  dataIndex: keyof T,
  options?: { format?: 'full' | 'relative'; width?: number; sorter?: boolean }
): ColumnType<T> => ({
  title,
  dataIndex: dataIndex as string,
  width: options?.width || 180,
  sorter: options?.sorter,
  render: (date: string) => {
    if (!date) return '-'
    return options?.format === 'relative'
      ? formatDateRelative(date)
      : formatDateTime(date)
  },
})

export const createStatusColumn = <T,>(
  title: string,
  dataIndex: keyof T,
  statusConfig: Record<string, { color: string; label?: string }>,
  options?: { width?: number }
): ColumnType<T> => ({
  title,
  dataIndex: dataIndex as string,
  width: options?.width || 120,
  render: (status: string) => {
    const config = statusConfig[status] || { color: 'default', label: status }
    return <Tag color={config.color}>{config.label || status}</Tag>
  },
})

export const createFileSizeColumn = <T,>(
  title: string,
  dataIndex: keyof T,
  options?: { width?: number; sorter?: boolean }
): ColumnType<T> => ({
  title,
  dataIndex: dataIndex as string,
  width: options?.width || 120,
  sorter: options?.sorter,
  render: (bytes: number) => (bytes ? formatFileSize(bytes) : '-'),
})
```

**Impact:** Eliminates ~50 lines per list page

---

#### Task 3.3: Consolidate Empty State Components
**Actions:**
1. Remove `frontend/src/components/common/EmptyState.tsx`
2. Keep `frontend/src/components/common/EnhancedEmptyState.tsx`
3. Rename to `EmptyState.tsx` (remove "Enhanced" prefix)
4. Update all imports

**Impact:** Single empty state component

---

### Phase 4: Store Optimization (Priority 4)

#### Task 4.1: Create Generic Store Factory
**New file:** `frontend/src/store/crud.factory.ts`

**Content:**
```typescript
import { create } from 'zustand'
import { message } from 'antd'
import { CrudService } from '@/services/crud.factory'
import { ListParams } from '@/types/common'

export interface CrudStoreState<T> {
  items: T[]
  current: T | null
  loading: boolean
  error: string | null
}

export interface CrudStoreActions<T, CreateT = Partial<T>, UpdateT = Partial<T>> {
  fetchItems: (params?: ListParams) => Promise<void>
  fetchById: (id: string) => Promise<void>
  createItem: (data: CreateT) => Promise<void>
  updateItem: (id: string, data: UpdateT) => Promise<void>
  deleteItem: (id: string) => Promise<void>
  setCurrent: (item: T | null) => void
  reset: () => void
}

export type CrudStore<T, CreateT = Partial<T>, UpdateT = Partial<T>> =
  CrudStoreState<T> & CrudStoreActions<T, CreateT, UpdateT>

/**
 * Create a CRUD store for a given service
 */
export function createCrudStore<T extends { id: string }, CreateT = Partial<T>, UpdateT = Partial<T>>(
  service: CrudService<T, CreateT, UpdateT>,
  entityName: string
) {
  return create<CrudStore<T, CreateT, UpdateT>>((set, get) => ({
    // State
    items: [],
    current: null,
    loading: false,
    error: null,

    // Actions
    fetchItems: async (params?: ListParams) => {
      set({ loading: true, error: null })
      try {
        const response = await service.getAll(params)
        set({ items: response.items, loading: false })
      } catch (error: any) {
        set({ error: error.message, loading: false })
        message.error(`Failed to fetch ${entityName}: ${error.message}`)
      }
    },

    fetchById: async (id: string) => {
      set({ loading: true, error: null })
      try {
        const item = await service.getById(id)
        set({ current: item, loading: false })
      } catch (error: any) {
        set({ error: error.message, loading: false })
        message.error(`Failed to fetch ${entityName}: ${error.message}`)
      }
    },

    createItem: async (data: CreateT) => {
      set({ loading: true, error: null })
      try {
        const newItem = await service.create(data)
        set((state) => ({
          items: [newItem, ...state.items],
          loading: false,
        }))
        message.success(`${entityName} created successfully`)
      } catch (error: any) {
        set({ error: error.message, loading: false })
        message.error(`Failed to create ${entityName}: ${error.message}`)
        throw error
      }
    },

    updateItem: async (id: string, data: UpdateT) => {
      set({ loading: true, error: null })
      try {
        const updatedItem = await service.update(id, data)
        set((state) => ({
          items: state.items.map((item) => (item.id === id ? updatedItem : item)),
          current: state.current?.id === id ? updatedItem : state.current,
          loading: false,
        }))
        message.success(`${entityName} updated successfully`)
      } catch (error: any) {
        set({ error: error.message, loading: false })
        message.error(`Failed to update ${entityName}: ${error.message}`)
        throw error
      }
    },

    deleteItem: async (id: string) => {
      set({ loading: true, error: null })
      try {
        await service.delete(id)
        set((state) => ({
          items: state.items.filter((item) => item.id !== id),
          current: state.current?.id === id ? null : state.current,
          loading: false,
        }))
        message.success(`${entityName} deleted successfully`)
      } catch (error: any) {
        set({ error: error.message, loading: false })
        message.error(`Failed to delete ${entityName}: ${error.message}`)
        throw error
      }
    },

    setCurrent: (item: T | null) => set({ current: item }),

    reset: () => set({ items: [], current: null, loading: false, error: null }),
  }))
}
```

**Impact:** Reduces store code from ~150 lines to ~10 lines each

---

## Testing Strategy

After each phase of refactoring:

1. **Unit Tests**: Test new utilities and factories
2. **Integration Tests**: Test refactored services and stores
3. **Manual Testing**: Verify UI functionality in each list page
4. **Regression Testing**: Ensure no features broken

---

## Success Metrics

- **Code Reduction**: ~2,000+ lines removed
- **File Count**: ~5 files removed, ~10 new utility files created
- **Maintainability**: 80% reduction in duplicate code
- **Type Safety**: 100% TypeScript coverage maintained
- **Test Coverage**: Maintain 60%+ frontend coverage
- **Build Size**: Reduce bundle size by ~5-10%

---

## Risk Mitigation

1. **Incremental Changes**: Refactor one module at a time
2. **Git Commits**: Commit after each successful refactor
3. **Backwards Compatibility**: Maintain old patterns during transition
4. **Testing**: Test each refactor before moving to next
5. **Rollback Plan**: Git history allows easy rollback

---

## Timeline

- **Phase 1**: 2-3 days (Critical infrastructure)
- **Phase 2**: 2-3 days (Service layer)
- **Phase 3**: 3-4 days (Components)
- **Phase 4**: 2-3 days (Stores)
- **Testing**: 1-2 days (Final verification)

**Total**: 10-15 days

---

## Conclusion

This cleanup plan addresses all identified redundancies systematically, prioritizing high-impact changes first. The result will be a more maintainable, consistent, and efficient codebase ready for production deployment.
