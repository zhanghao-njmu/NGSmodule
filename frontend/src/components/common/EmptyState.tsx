/**
 * Empty State Component - Consistent empty state display
 */
import React from 'react'
import { Empty, Button } from 'antd'
import { InboxOutlined } from '@ant-design/icons'

interface EmptyStateProps {
  icon?: React.ReactNode
  title?: string
  description?: string
  actionText?: string
  onAction?: () => void
  style?: React.CSSProperties
}

export const EmptyState: React.FC<EmptyStateProps> = ({
  icon,
  title = 'No Data',
  description = 'There are no records to display',
  actionText,
  onAction,
  style,
}) => {
  return (
    <div className="empty-state" style={style}>
      <Empty
        image={
          icon || (
            <InboxOutlined
              style={{
                fontSize: 64,
                color: 'var(--text-tertiary)',
              }}
            />
          )
        }
        imageStyle={{
          height: 80,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
        }}
        description={
          <div>
            <div className="empty-state-title">{title}</div>
            <div className="empty-state-description">{description}</div>
          </div>
        }
      >
        {actionText && onAction && (
          <Button type="primary" onClick={onAction}>
            {actionText}
          </Button>
        )}
      </Empty>
    </div>
  )
}

export default EmptyState
