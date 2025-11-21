/**
 * Page Header Component - Consistent header with filters and actions
 */
import React from 'react'
import { Card, Space } from 'antd'

interface PageHeaderProps {
  left?: React.ReactNode
  right?: React.ReactNode
  style?: React.CSSProperties
  className?: string
}

export const PageHeader: React.FC<PageHeaderProps> = ({
  left,
  right,
  style,
  className,
}) => {
  return (
    <Card style={{ marginBottom: 16, ...style }} className={className}>
      <Space style={{ width: '100%', justifyContent: 'space-between', flexWrap: 'wrap' }}>
        {left && <Space wrap>{left}</Space>}
        {right && <Space wrap>{right}</Space>}
      </Space>
    </Card>
  )
}

export default PageHeader
