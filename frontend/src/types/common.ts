/**
 * Common Types
 * Shared type definitions used across the application
 */

/**
 * Generic paginated response structure
 * Used for list endpoints that return paginated data
 *
 * @example
 * const response: PaginatedResponse<Project> = await api.get('/projects')
 * console.log(response.total, response.items)
 */
export interface PaginatedResponse<T> {
  /** Total number of items available */
  total: number
  /** Array of items for current page */
  items: T[]
  /** Current page number (optional, 1-indexed) */
  page?: number
  /** Number of items per page (optional) */
  page_size?: number
  /** Offset for pagination (optional, 0-indexed) */
  skip?: number
  /** Maximum number of items to return (optional) */
  limit?: number
}

/**
 * Generic list/query parameters
 * Used when fetching paginated lists
 *
 * @example
 * const params: ListParams = { page: 1, page_size: 20, search: 'RNA-seq' }
 * const projects = await projectService.getAll(params)
 */
export interface ListParams {
  /** Page number (1-indexed) */
  page?: number
  /** Number of items per page */
  page_size?: number
  /** Offset (0-indexed, alternative to page) */
  skip?: number
  /** Maximum number of items (alternative to page_size) */
  limit?: number
  /** Search query string */
  search?: string
  /** Field to sort by */
  sort_by?: string
  /** Sort direction */
  sort_order?: 'asc' | 'desc' | 'ascend' | 'descend'
  /** Additional filters */
  [key: string]: any
}

/**
 * Generic API success response
 * Used for endpoints that return a single data object with metadata
 */
export interface ApiResponse<T> {
  /** Indicates if the request was successful */
  success: boolean
  /** The response data */
  data: T
  /** Optional message (success/info message) */
  message?: string
  /** Optional metadata */
  meta?: Record<string, any>
}

/**
 * Generic API error response
 * Standardized error structure from backend
 */
export interface ApiError {
  /** Error message */
  message: string
  /** Error code (e.g., 'VALIDATION_ERROR', 'NOT_FOUND') */
  code?: string
  /** HTTP status code */
  status?: number
  /** Additional error details */
  details?: Record<string, any>
  /** Field-specific errors (for validation errors) */
  errors?: Array<{
    field: string
    message: string
    code?: string
  }>
}

/**
 * Generic statistics structure
 * Used for dashboard and analytics data
 */
export interface Stats {
  /** Total count */
  total: number
  /** Active/in-progress count */
  active?: number
  /** Completed count */
  completed?: number
  /** Failed/error count */
  failed?: number
  /** Pending count */
  pending?: number
  /** Success rate (0-1) */
  success_rate?: number
  /** Additional metrics */
  [key: string]: number | undefined
}

/**
 * Generic filter configuration
 * Used in filter components
 */
export interface FilterConfig {
  /** Filter field name */
  field: string
  /** Display label */
  label: string
  /** Filter type */
  type: 'select' | 'date' | 'dateRange' | 'number' | 'text' | 'checkbox'
  /** Options for select type */
  options?: Array<{ label: string; value: any }>
  /** Placeholder text */
  placeholder?: string
  /** Default value */
  defaultValue?: any
}

/**
 * Generic sort configuration
 * Used in sortable tables
 */
export interface SortConfig {
  /** Field being sorted */
  field: string
  /** Sort direction */
  order: 'asc' | 'desc' | 'ascend' | 'descend'
}

/**
 * Generic ID-based entity
 * Base interface for all entities with an ID
 */
export interface Entity {
  /** Unique identifier */
  id: string
  /** Creation timestamp */
  created_at?: string
  /** Last update timestamp */
  updated_at?: string
}

/**
 * Generic action result
 * Used for operations that don't return specific data
 */
export interface ActionResult {
  /** Whether the action succeeded */
  success: boolean
  /** Result message */
  message: string
  /** Optional result data */
  data?: any
}

/**
 * Generic option type
 * Used in selects, dropdowns, etc.
 */
export interface Option<T = any> {
  /** Display label */
  label: string
  /** Option value */
  value: T
  /** Whether option is disabled */
  disabled?: boolean
  /** Optional icon */
  icon?: React.ReactNode
  /** Additional metadata */
  meta?: Record<string, any>
}

/**
 * File upload info
 * Used for file upload operations
 */
export interface UploadInfo {
  /** File name */
  name: string
  /** File size in bytes */
  size: number
  /** MIME type */
  type: string
  /** Upload status */
  status: 'uploading' | 'done' | 'error'
  /** Upload progress (0-100) */
  percent?: number
  /** Server response after upload */
  response?: any
  /** Error message if failed */
  error?: string
}

/**
 * Generic key-value pair
 */
export interface KeyValue<K = string, V = any> {
  key: K
  value: V
}

/**
 * Generic tree node structure
 * Used for tree/hierarchical data
 */
export interface TreeNode<T = any> {
  /** Node key */
  key: string
  /** Node title/label */
  title: string
  /** Child nodes */
  children?: TreeNode<T>[]
  /** Whether node is a leaf */
  isLeaf?: boolean
  /** Whether node is disabled */
  disabled?: boolean
  /** Additional node data */
  data?: T
}

/**
 * Breadcrumb item
 */
export interface BreadcrumbItem {
  /** Display title */
  title: string
  /** Navigation path */
  path?: string
  /** Icon component */
  icon?: React.ReactNode
}

/**
 * Tab item configuration
 */
export interface TabItem {
  /** Tab key */
  key: string
  /** Tab label */
  label: string
  /** Tab content */
  children: React.ReactNode
  /** Whether tab is disabled */
  disabled?: boolean
  /** Icon */
  icon?: React.ReactNode
  /** Close button */
  closable?: boolean
}

/**
 * Menu item configuration
 */
export interface MenuItem {
  /** Item key */
  key: string
  /** Display label */
  label: string
  /** Icon */
  icon?: React.ReactNode
  /** Navigation path */
  path?: string
  /** Sub-menu items */
  children?: MenuItem[]
  /** Whether item is disabled */
  disabled?: boolean
  /** Required permissions */
  permissions?: string[]
}

/**
 * Notification item
 */
export interface NotificationItem {
  /** Notification ID */
  id: string
  /** Notification type */
  type: 'info' | 'success' | 'warning' | 'error'
  /** Title */
  title: string
  /** Message content */
  message: string
  /** Timestamp */
  timestamp: string
  /** Whether notification is read */
  read: boolean
  /** Action URL */
  actionUrl?: string
  /** Action label */
  actionLabel?: string
}
