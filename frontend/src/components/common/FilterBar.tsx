/**
 * FilterBar - 通用过滤栏组件
 * 支持搜索、下拉选择、日期范围等常见过滤场景
 */
import React from 'react'
import { Input, Select, DatePicker, Space, Button } from 'antd'
import { SearchOutlined, ReloadOutlined } from '@ant-design/icons'
import type { SelectProps } from 'antd'

const { Search } = Input
const { RangePicker } = DatePicker

export interface FilterConfig {
  type: 'search' | 'select' | 'dateRange'
  key: string
  placeholder?: string
  label?: string
  options?: SelectProps['options']
  allowClear?: boolean
  showSearch?: boolean
  width?: number | string
}

export interface FilterBarProps {
  filters: FilterConfig[]
  onFilterChange?: (key: string, value: any) => void
  onReset?: () => void
  showReset?: boolean
  values?: Record<string, any>
}

export const FilterBar: React.FC<FilterBarProps> = ({
  filters,
  onFilterChange,
  onReset,
  showReset = true,
  values = {},
}) => {
  const handleFilterChange = (key: string, value: any) => {
    onFilterChange?.(key, value)
  }

  const handleReset = () => {
    onReset?.()
  }

  const renderFilter = (filter: FilterConfig) => {
    const commonProps = {
      key: filter.key,
      placeholder: filter.placeholder,
      allowClear: filter.allowClear ?? true,
      style: { width: filter.width || 200 },
      value: values[filter.key],
    }

    switch (filter.type) {
      case 'search':
        return (
          <Search
            {...commonProps}
            prefix={<SearchOutlined />}
            onSearch={(value) => handleFilterChange(filter.key, value)}
            onChange={(e) => handleFilterChange(filter.key, e.target.value)}
            enterButton={false}
          />
        )

      case 'select':
        return (
          <Select
            {...commonProps}
            options={filter.options}
            showSearch={filter.showSearch}
            onChange={(value) => handleFilterChange(filter.key, value)}
            filterOption={(input, option) =>
              (option?.label?.toString() ?? '')
                .toLowerCase()
                .includes(input.toLowerCase())
            }
          />
        )

      case 'dateRange':
        return (
          <RangePicker
            {...commonProps}
            onChange={(dates) => handleFilterChange(filter.key, dates)}
          />
        )

      default:
        return null
    }
  }

  return (
    <Space size="middle" wrap style={{ marginBottom: 16 }}>
      {filters.map((filter) => (
        <div key={filter.key}>
          {filter.label && (
            <span style={{ marginRight: 8, color: 'var(--text-secondary)' }}>
              {filter.label}:
            </span>
          )}
          {renderFilter(filter)}
        </div>
      ))}
      {showReset && (
        <Button icon={<ReloadOutlined />} onClick={handleReset}>
          Reset
        </Button>
      )}
    </Space>
  )
}

export default FilterBar
