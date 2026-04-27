/**
 * TableSkeleton - 表格骨架屏
 * 用于表格加载时的占位显示
 */
import type React from 'react'
import { Skeleton, Table } from 'antd'
import type { ColumnsType } from 'antd/es/table'

export interface TableSkeletonProps {
  columns?: number
  rows?: number
  avatar?: boolean
}

export const TableSkeleton: React.FC<TableSkeletonProps> = ({ columns = 5, rows = 8, avatar = false }) => {
  // 生成骨架屏列
  const skeletonColumns: ColumnsType<any> = Array.from({ length: columns }, (_, index) => ({
    key: `col-${index}`,
    title: <Skeleton.Input active size="small" style={{ width: 80 }} />,
    dataIndex: `col-${index}`,
    render: () => <Skeleton active paragraph={false} title={{ width: '100%' }} avatar={avatar && index === 0} />,
  }))

  // 生成骨架屏数据
  const skeletonData = Array.from({ length: rows }, (_, index) => ({
    key: `row-${index}`,
  }))

  return <Table columns={skeletonColumns} dataSource={skeletonData} pagination={false} showHeader={true} />
}

export default TableSkeleton
