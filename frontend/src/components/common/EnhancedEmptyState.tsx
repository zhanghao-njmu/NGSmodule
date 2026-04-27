/**
 * EnhancedEmptyState - 增强版空状态组件
 * 支持多种场景和自定义操作
 */
import React from 'react'
import { Empty, Button, Typography } from 'antd'
import {
  FileOutlined,
  InboxOutlined,
  FolderOpenOutlined,
  SearchOutlined,
  LockOutlined,
  ExclamationCircleOutlined,
} from '@ant-design/icons'
import type { EmptyProps } from 'antd'

const { Text, Paragraph } = Typography

export type EmptyStateType = 'noData' | 'noSearchResults' | 'noPermission' | 'error' | 'empty' | 'custom'

export interface EnhancedEmptyStateProps extends Omit<EmptyProps, 'image'> {
  type?: EmptyStateType
  title?: string
  description?: string | React.ReactNode
  action?: {
    text: string
    onClick: () => void
    icon?: React.ReactNode
    type?: 'primary' | 'default' | 'dashed'
  }
  customIcon?: React.ReactNode
  size?: 'small' | 'default' | 'large'
}

const iconMap: Record<EmptyStateType, React.ReactNode> = {
  noData: <InboxOutlined style={{ fontSize: 64, color: 'var(--color-gray-300)' }} />,
  noSearchResults: <SearchOutlined style={{ fontSize: 64, color: 'var(--color-gray-300)' }} />,
  noPermission: <LockOutlined style={{ fontSize: 64, color: 'var(--color-warning)' }} />,
  error: <ExclamationCircleOutlined style={{ fontSize: 64, color: 'var(--color-error)' }} />,
  empty: <FolderOpenOutlined style={{ fontSize: 64, color: 'var(--color-gray-300)' }} />,
  custom: <FileOutlined style={{ fontSize: 64, color: 'var(--color-gray-300)' }} />,
}

const titleMap: Record<EmptyStateType, string> = {
  noData: 'No Data',
  noSearchResults: 'No Results Found',
  noPermission: 'Access Denied',
  error: 'Something Went Wrong',
  empty: 'No Items',
  custom: 'Empty',
}

const descriptionMap: Record<EmptyStateType, string> = {
  noData: 'There is no data to display yet.',
  noSearchResults: 'Try adjusting your search criteria.',
  noPermission: 'You do not have permission to view this content.',
  error: 'An error occurred while loading data.',
  empty: 'This folder is empty.',
  custom: '',
}

export const EnhancedEmptyState: React.FC<EnhancedEmptyStateProps> = ({
  type = 'noData',
  title,
  description,
  action,
  customIcon,
  size = 'default',
  ...rest
}) => {
  const imageSize = size === 'small' ? 48 : size === 'large' ? 80 : 64
  const icon = customIcon || iconMap[type]

  // Clone icon with size
  const sizedIcon = React.isValidElement(icon)
    ? React.cloneElement(icon as React.ReactElement<any>, {
        style: { ...icon.props.style, fontSize: imageSize },
      })
    : icon

  const finalTitle = title || titleMap[type]
  const finalDescription = description || descriptionMap[type]

  return (
    <Empty
      image={sizedIcon}
      description={
        <div style={{ marginTop: size === 'small' ? 8 : 16 }}>
          {finalTitle && (
            <Text
              strong
              style={{
                display: 'block',
                fontSize: size === 'small' ? 14 : 16,
                color: 'var(--text-primary)',
                marginBottom: 4,
              }}
            >
              {finalTitle}
            </Text>
          )}
          {finalDescription && (
            <Paragraph
              type="secondary"
              style={{
                margin: 0,
                fontSize: size === 'small' ? 12 : 14,
              }}
            >
              {finalDescription}
            </Paragraph>
          )}
        </div>
      }
      {...rest}
    >
      {action && (
        <Button
          type={action.type || 'primary'}
          icon={action.icon}
          onClick={action.onClick}
          style={{ marginTop: size === 'small' ? 12 : 16 }}
        >
          {action.text}
        </Button>
      )}
    </Empty>
  )
}

export default EnhancedEmptyState
