/**
 * PageSkeleton - 页面级骨架屏
 * 用于页面加载时的占位显示
 */
import React from 'react'
import { Card, Skeleton, Space } from 'antd'

export interface PageSkeletonProps {
  hasHeader?: boolean
  hasFilters?: boolean
  rows?: number
  columns?: number
}

export const PageSkeleton: React.FC<PageSkeletonProps> = ({
  hasHeader = true,
  hasFilters = true,
  rows = 5,
  columns = 4,
}) => {
  return (
    <div style={{ padding: '24px' }}>
      {/* Header Skeleton */}
      {hasHeader && (
        <div
          style={{
            display: 'flex',
            justifyContent: 'space-between',
            marginBottom: 24,
          }}
        >
          <Skeleton.Input active style={{ width: 200, height: 36 }} />
          <Skeleton.Button active style={{ width: 120, height: 36 }} />
        </div>
      )}

      {/* Filters Skeleton */}
      {hasFilters && (
        <Space size="middle" style={{ marginBottom: 16 }}>
          <Skeleton.Input active style={{ width: 200, height: 36 }} />
          <Skeleton.Input active style={{ width: 150, height: 36 }} />
          <Skeleton.Button active style={{ width: 80, height: 36 }} />
        </Space>
      )}

      {/* Table Skeleton */}
      <Card>
        <Skeleton active paragraph={{ rows }} />
      </Card>
    </div>
  )
}

export default PageSkeleton
