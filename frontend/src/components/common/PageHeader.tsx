/**
 * Page Header Component - Consistent header with filters and actions
 */
import type React from 'react'
import { Card, Space, Typography } from 'antd'

const { Title, Text } = Typography

export interface PageHeaderProps {
  title?: string | React.ReactNode
  subtitle?: string | React.ReactNode
  icon?: React.ReactNode
  extra?: React.ReactNode
  left?: React.ReactNode
  right?: React.ReactNode
  style?: React.CSSProperties
  className?: string
}

export const PageHeader: React.FC<PageHeaderProps> = ({
  title,
  subtitle,
  icon,
  extra,
  left,
  right,
  style,
  className,
}) => {
  return (
    <Card style={{ marginBottom: 16, ...style }} className={className}>
      <Space style={{ width: '100%', justifyContent: 'space-between', flexWrap: 'wrap' }}>
        {(title || subtitle || icon || left) && (
          <Space wrap>
            {icon}
            <div>
              {title && (
                <Title level={4} style={{ margin: 0 }}>
                  {title}
                </Title>
              )}
              {subtitle && <Text type="secondary">{subtitle}</Text>}
            </div>
            {left}
          </Space>
        )}
        {(right || extra) && <Space wrap>{right || extra}</Space>}
      </Space>
    </Card>
  )
}

export default PageHeader
