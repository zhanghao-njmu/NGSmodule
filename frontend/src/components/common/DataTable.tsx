/**
 * Data Table Component - Consistent table with pagination
 */
import type React from 'react'
import { Card, Table, Empty } from 'antd'
import type { TableProps } from 'antd'
import { FileTextOutlined } from '@ant-design/icons'

interface DataTableProps<T> extends Omit<TableProps<T>, 'title'> {
  title?: React.ReactNode
  extra?: React.ReactNode
  cardStyle?: React.CSSProperties
  emptyText?: string
  emptyDescription?: string
}

export function DataTable<T extends object>({
  title,
  extra,
  cardStyle,
  emptyText = 'No Data',
  emptyDescription = 'There are no records to display',
  pagination = {
    pageSize: 20,
    showSizeChanger: true,
    showTotal: (total) => `Total ${total} records`,
  },
  locale = {
    emptyText: (
      <Empty
        image={<FileTextOutlined style={{ fontSize: 64, color: 'var(--text-tertiary)' }} />}
        description={
          <div>
            <div style={{ fontWeight: 500, marginBottom: 4 }}>{emptyText}</div>
            <div style={{ color: 'var(--text-tertiary)', fontSize: 14 }}>{emptyDescription}</div>
          </div>
        }
      />
    ),
  },
  ...tableProps
}: DataTableProps<T>) {
  return (
    <Card title={title} extra={extra} style={cardStyle}>
      <Table {...tableProps} pagination={pagination} locale={locale} className="data-table" />
    </Card>
  )
}

export default DataTable
