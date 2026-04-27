/**
 * Status Tag Component - Consistent status display with icons
 */
import type React from 'react'
import { Tag } from 'antd'
import type { TagProps } from 'antd'
import {
  CheckCircleOutlined,
  CloseCircleOutlined,
  ClockCircleOutlined,
  SyncOutlined,
  StopOutlined,
  InboxOutlined,
} from '@ant-design/icons'

export type StatusType =
  | 'success'
  | 'error'
  | 'warning'
  | 'processing'
  | 'default'
  | 'pending'
  | 'running'
  | 'completed'
  | 'failed'
  | 'cancelled'
  | 'active'
  | 'archived'
  | 'inactive'

interface StatusConfig {
  color: string
  icon?: React.ReactNode
}

const statusConfigMap: Record<string, StatusConfig> = {
  // Task statuses
  pending: { color: 'default', icon: <ClockCircleOutlined /> },
  running: { color: 'processing', icon: <SyncOutlined spin /> },
  completed: { color: 'success', icon: <CheckCircleOutlined /> },
  failed: { color: 'error', icon: <CloseCircleOutlined /> },
  cancelled: { color: 'warning', icon: <StopOutlined /> },

  // Project statuses
  active: { color: 'success', icon: <CheckCircleOutlined /> },
  archived: { color: 'default', icon: <InboxOutlined /> },
  inactive: { color: 'default', icon: <StopOutlined /> },

  // Generic statuses
  success: { color: 'success', icon: <CheckCircleOutlined /> },
  error: { color: 'error', icon: <CloseCircleOutlined /> },
  warning: { color: 'warning', icon: <ClockCircleOutlined /> },
  processing: { color: 'processing', icon: <SyncOutlined spin /> },
  default: { color: 'default' },
}

interface StatusTagProps extends Omit<TagProps, 'color' | 'icon'> {
  status: StatusType | string
  showIcon?: boolean
  text?: string
}

export const StatusTag: React.FC<StatusTagProps> = ({ status, showIcon = true, text, children, ...tagProps }) => {
  const config = statusConfigMap[status.toLowerCase()] || statusConfigMap.default
  const displayText = text || children || status.toUpperCase()

  return (
    <Tag color={config.color} icon={showIcon ? config.icon : undefined} {...tagProps}>
      {displayText}
    </Tag>
  )
}

export default StatusTag
