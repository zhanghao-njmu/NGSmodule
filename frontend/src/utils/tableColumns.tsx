/**
 * Table Column Factories
 * Reusable column definitions for Ant Design tables
 * Eliminates code duplication across list pages
 */

import React from 'react'
import { Space, Button, Tooltip, Popconfirm, Tag, Typography } from 'antd'
import {
  EditOutlined,
  DeleteOutlined,
  EyeOutlined,
  CopyOutlined,
  DownloadOutlined,
  PlayCircleOutlined,
  StopOutlined,
} from '@ant-design/icons'
import type { ColumnType } from 'antd/es/table'
import { formatDateTime, formatDateRelative, formatFileSize, formatPercent, truncate } from './format'

const { Text, Link } = Typography

/**
 * Status color configuration
 */
export interface StatusConfig {
  color: string
  label?: string
}

/**
 * Action button configuration
 */
export interface ActionButton<T> {
  icon: React.ReactNode
  tooltip: string
  onClick: (record: T) => void
  danger?: boolean
  disabled?: (record: T) => boolean
  visible?: (record: T) => boolean
  confirm?: {
    title: string
    okText?: string
    cancelText?: string
  }
}

/**
 * Create an actions column with edit/delete buttons
 *
 * @param onEdit - Edit button click handler
 * @param onDelete - Delete button click handler
 * @param onView - Optional view button click handler
 * @param options - Additional options
 * @returns Table column configuration
 *
 * @example
 * createActionColumn(
 *   (record) => showEditModal(record),
 *   (record) => handleDelete(record.id),
 *   (record) => navigate(`/projects/${record.id}`)
 * )
 */
export function createActionColumn<T extends { id: string }>(
  onEdit?: (record: T) => void,
  onDelete?: (record: T) => void,
  onView?: (record: T) => void,
  options?: {
    width?: number
    fixed?: 'left' | 'right'
    deleteConfirmTitle?: string
    deleteOkText?: string
    deleteCancelText?: string
    showEdit?: boolean
    showDelete?: boolean
    showView?: boolean
  }
): ColumnType<T> {
  const {
    width = 120,
    fixed = 'right',
    deleteConfirmTitle = 'Are you sure you want to delete this item?',
    deleteOkText = 'Yes',
    deleteCancelText = 'No',
    showEdit = true,
    showDelete = true,
    showView = false,
  } = options || {}

  return {
    title: 'Actions',
    key: 'actions',
    width,
    fixed,
    render: (_, record) => (
      <Space size="small">
        {showView && onView && (
          <Tooltip title="View">
            <Button
              type="text"
              size="small"
              icon={<EyeOutlined />}
              onClick={() => onView(record)}
            />
          </Tooltip>
        )}

        {showEdit && onEdit && (
          <Tooltip title="Edit">
            <Button
              type="text"
              size="small"
              icon={<EditOutlined />}
              onClick={() => onEdit(record)}
            />
          </Tooltip>
        )}

        {showDelete && onDelete && (
          <Popconfirm
            title={deleteConfirmTitle}
            onConfirm={() => onDelete(record)}
            okText={deleteOkText}
            cancelText={deleteCancelText}
          >
            <Tooltip title="Delete">
              <Button
                type="text"
                size="small"
                danger
                icon={<DeleteOutlined />}
              />
            </Tooltip>
          </Popconfirm>
        )}
      </Space>
    ),
  }
}

/**
 * Create a custom actions column with flexible button configuration
 *
 * @param buttons - Array of action button configurations
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createCustomActionColumn([
 *   {
 *     icon: <PlayCircleOutlined />,
 *     tooltip: 'Run',
 *     onClick: (record) => runPipeline(record.id),
 *   },
 *   {
 *     icon: <DownloadOutlined />,
 *     tooltip: 'Download',
 *     onClick: (record) => downloadResults(record.id),
 *   }
 * ])
 */
export function createCustomActionColumn<T>(
  buttons: ActionButton<T>[],
  options?: { width?: number; fixed?: 'left' | 'right' }
): ColumnType<T> {
  const { width = 150, fixed = 'right' } = options || {}

  return {
    title: 'Actions',
    key: 'actions',
    width,
    fixed,
    render: (_, record) => (
      <Space size="small">
        {buttons.map((button, index) => {
          // Check if button should be visible
          if (button.visible && !button.visible(record)) {
            return null
          }

          // Check if button should be disabled
          const isDisabled = button.disabled ? button.disabled(record) : false

          const buttonElement = (
            <Tooltip title={button.tooltip} key={index}>
              <Button
                type="text"
                size="small"
                icon={button.icon}
                onClick={() => button.onClick(record)}
                danger={button.danger}
                disabled={isDisabled}
              />
            </Tooltip>
          )

          // Wrap with Popconfirm if confirm config provided
          if (button.confirm) {
            return (
              <Popconfirm
                key={index}
                title={button.confirm.title}
                onConfirm={() => button.onClick(record)}
                okText={button.confirm.okText}
                cancelText={button.confirm.cancelText}
              >
                {buttonElement}
              </Popconfirm>
            )
          }

          return buttonElement
        })}
      </Space>
    ),
  }
}

/**
 * Create a date column with formatting
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createDateColumn('Created At', 'created_at')
 * createDateColumn('Updated', 'updated_at', { format: 'relative' })
 */
export function createDateColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    format?: 'full' | 'relative' | 'custom'
    customFormat?: string
    width?: number
    sorter?: boolean
    defaultSortOrder?: 'ascend' | 'descend'
  }
): ColumnType<T> {
  const {
    format = 'full',
    customFormat,
    width = 180,
    sorter = false,
    defaultSortOrder,
  } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    sorter: sorter
      ? (a: any, b: any) => {
          const dateA = new Date(a[dataIndex]).getTime()
          const dateB = new Date(b[dataIndex]).getTime()
          return dateA - dateB
        }
      : undefined,
    defaultSortOrder,
    render: (date: string) => {
      if (!date) return <Text type="secondary">-</Text>

      if (format === 'relative') {
        return formatDateRelative(date)
      }

      if (format === 'custom' && customFormat) {
        return formatDateTime(date, customFormat)
      }

      return formatDateTime(date)
    },
  }
}

/**
 * Create a status column with colored tags
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param statusConfig - Status color configuration
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createStatusColumn('Status', 'status', {
 *   active: { color: 'green', label: 'Active' },
 *   inactive: { color: 'gray', label: 'Inactive' },
 *   pending: { color: 'orange' }
 * })
 */
export function createStatusColumn<T>(
  title: string,
  dataIndex: keyof T,
  statusConfig: Record<string, StatusConfig>,
  options?: {
    width?: number
    filters?: boolean
  }
): ColumnType<T> {
  const { width = 120, filters = false } = options || {}

  // Generate filters from status config if enabled
  const columnFilters = filters
    ? Object.entries(statusConfig).map(([status, config]) => ({
        text: config.label || status,
        value: status,
      }))
    : undefined

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    filters: columnFilters,
    onFilter: filters
      ? (value: any, record: any) => record[dataIndex] === value
      : undefined,
    render: (status: string) => {
      const config = statusConfig[status] || { color: 'default', label: status }
      return (
        <Tag color={config.color}>
          {config.label || status?.toUpperCase() || 'UNKNOWN'}
        </Tag>
      )
    },
  }
}

/**
 * Create a file size column
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createFileSizeColumn('Size', 'file_size')
 */
export function createFileSizeColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    sorter?: boolean
  }
): ColumnType<T> {
  const { width = 120, sorter = false } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    sorter: sorter
      ? (a: any, b: any) => a[dataIndex] - b[dataIndex]
      : undefined,
    render: (bytes: number) =>
      bytes ? formatFileSize(bytes) : <Text type="secondary">-</Text>,
  }
}

/**
 * Create a percentage column
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createPercentColumn('Progress', 'progress')
 */
export function createPercentColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    decimals?: number
    sorter?: boolean
    colorize?: boolean // Color based on value (green > 80%, yellow > 60%, red < 60%)
  }
): ColumnType<T> {
  const { width = 100, decimals = 1, sorter = false, colorize = false } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    sorter: sorter
      ? (a: any, b: any) => a[dataIndex] - b[dataIndex]
      : undefined,
    render: (value: number) => {
      if (value === null || value === undefined) {
        return <Text type="secondary">-</Text>
      }

      const formatted = formatPercent(value, decimals)

      if (!colorize) {
        return formatted
      }

      // Colorize based on value
      const percentage = value * 100
      let color: string

      if (percentage >= 80) {
        color = '#52c41a' // Green
      } else if (percentage >= 60) {
        color = '#faad14' // Yellow
      } else {
        color = '#ff4d4f' // Red
      }

      return <Text style={{ color }}>{formatted}</Text>
    },
  }
}

/**
 * Create a link column
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param onClick - Link click handler
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createLinkColumn(
 *   'Project Name',
 *   'name',
 *   (record) => navigate(`/projects/${record.id}`)
 * )
 */
export function createLinkColumn<T>(
  title: string,
  dataIndex: keyof T,
  onClick: (record: T) => void,
  options?: {
    width?: number
    ellipsis?: boolean
    maxLength?: number
  }
): ColumnType<T> {
  const { width, ellipsis = true, maxLength } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    ellipsis,
    render: (text: string, record: T) => {
      const displayText = maxLength ? truncate(text, maxLength) : text

      return (
        <Link onClick={() => onClick(record)}>
          {displayText || <Text type="secondary">-</Text>}
        </Link>
      )
    },
  }
}

/**
 * Create a text column with optional truncation
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createTextColumn('Description', 'description', { maxLength: 50 })
 */
export function createTextColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    ellipsis?: boolean
    maxLength?: number
    type?: 'default' | 'secondary' | 'success' | 'warning' | 'danger'
  }
): ColumnType<T> {
  const { width, ellipsis = false, maxLength, type = 'default' } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    ellipsis,
    render: (text: string) => {
      if (!text) return <Text type="secondary">-</Text>

      const displayText = maxLength ? truncate(text, maxLength) : text

      return <Text type={type}>{displayText}</Text>
    },
  }
}

/**
 * Create a boolean column (Yes/No or custom labels)
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createBooleanColumn('Active', 'is_active')
 * createBooleanColumn('Public', 'is_public', { labels: ['Public', 'Private'] })
 */
export function createBooleanColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    labels?: [string, string] // [true label, false label]
    colors?: [string, string] // [true color, false color]
  }
): ColumnType<T> {
  const {
    width = 100,
    labels = ['Yes', 'No'],
    colors = ['green', 'default'],
  } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    render: (value: boolean) => {
      if (value === null || value === undefined) {
        return <Text type="secondary">-</Text>
      }

      const label = value ? labels[0] : labels[1]
      const color = value ? colors[0] : colors[1]

      return <Tag color={color}>{label}</Tag>
    },
  }
}

/**
 * Create a number column with formatting
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createNumberColumn('Count', 'sample_count', { sorter: true })
 */
export function createNumberColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    decimals?: number
    sorter?: boolean
    prefix?: string
    suffix?: string
  }
): ColumnType<T> {
  const { width = 120, decimals = 0, sorter = false, prefix, suffix } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    sorter: sorter
      ? (a: any, b: any) => a[dataIndex] - b[dataIndex]
      : undefined,
    render: (value: number) => {
      if (value === null || value === undefined) {
        return <Text type="secondary">-</Text>
      }

      const formatted = value.toFixed(decimals)
      return `${prefix || ''}${formatted}${suffix || ''}`
    },
  }
}

/**
 * Create a tags column (array of strings as tags)
 *
 * @param title - Column title
 * @param dataIndex - Data field name
 * @param options - Column options
 * @returns Table column configuration
 *
 * @example
 * createTagsColumn('Labels', 'tags')
 */
export function createTagsColumn<T>(
  title: string,
  dataIndex: keyof T,
  options?: {
    width?: number
    color?: string
    maxTags?: number
  }
): ColumnType<T> {
  const { width = 200, color = 'blue', maxTags } = options || {}

  return {
    title,
    dataIndex: dataIndex as string,
    key: dataIndex as string,
    width,
    render: (tags: string[]) => {
      if (!tags || tags.length === 0) {
        return <Text type="secondary">-</Text>
      }

      const visibleTags = maxTags ? tags.slice(0, maxTags) : tags
      const remaining = maxTags && tags.length > maxTags ? tags.length - maxTags : 0

      return (
        <Space size={[0, 4]} wrap>
          {visibleTags.map((tag, index) => (
            <Tag key={index} color={color}>
              {tag}
            </Tag>
          ))}
          {remaining > 0 && <Tag>+{remaining} more</Tag>}
        </Space>
      )
    },
  }
}

export default {
  createActionColumn,
  createCustomActionColumn,
  createDateColumn,
  createStatusColumn,
  createFileSizeColumn,
  createPercentColumn,
  createLinkColumn,
  createTextColumn,
  createBooleanColumn,
  createNumberColumn,
  createTagsColumn,
}
