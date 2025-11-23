/**
 * Result List Page - Browse all analysis results
 */
import type React from 'react'
import { useEffect, useMemo, useState } from 'react'
import { useNavigate } from 'react-router-dom'
import {
  Card,
  Table,
  Tag,
  Button,
  Space,
  Tooltip,
  Typography,
  Row,
  Col,
  Statistic,
  Dropdown,
  type MenuProps,
} from 'antd'
import {
  EyeOutlined,
  BarChartOutlined,
  ExperimentOutlined,
  CheckCircleOutlined,
  WarningOutlined,
  CloseCircleOutlined,
  DownloadOutlined,
  MoreOutlined,
  ThunderboltOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import resultService from '@/services/result.service'
import { FilterBar, EnhancedEmptyState, PageSkeleton, FadeIn, StaggeredList, StatusTag } from '@/components/common'
import type { FilterConfig } from '@/components/common'
import { useAsync, useFilters, usePagination } from '@/hooks'
import { toast } from '@/utils/notification'
import type { Result } from '@/types/result'
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(relativeTime)

const { Title, Text } = Typography

export const ResultList: React.FC = () => {
  const navigate = useNavigate()
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([])

  // Filters using custom hook
  const { filters, setFilter, resetFilters } = useFilters({
    initialFilters: {
      search: '',
      result_type: 'all',
      task_id: undefined,
    },
  })

  // Pagination using custom hook
  const {
    pagination,
    onChange: handlePageChange,
    reset: resetPagination,
    skip,
    limit,
  } = usePagination({
    initialPageSize: 20,
  })

  // Load results with filters and pagination
  const {
    data: resultsData,
    loading,
    error,
    execute: loadResults,
    error,
  } = useAsync(
    async () => {
      const params: any = {
        skip,
        limit,
      }

      if (filters.result_type && filters.result_type !== 'all') {
        params.result_type = filters.result_type
      }

      if (filters.task_id) {
        params.task_id = filters.task_id
      }

      return resultService.getAll(params)
    },
    {
      immediate: false,
      onError: (err) => toast.error(`Failed to load results: ${err.message}`),
    },
  )

  // Load results when filters or pagination changes
  useEffect(() => {
    loadResults()
  }, [filters, pagination.current, pagination.pageSize])

  // Filter results by search term (client-side)
  const filteredResults = useMemo(() => {
    if (!resultsData?.results) {
      return []
    }

    if (!filters.search) {
      return resultsData.items
    }

    const searchLower = filters.search.toLowerCase()
    return resultsData.items.filter(
      (result: Result) =>
        result.result_type.toLowerCase().includes(searchLower) || result.id.toString().includes(searchLower),
    )
  }, [resultsData, filters.search])

  // Calculate statistics
  const stats = useMemo(() => {
    const results = resultsData?.results || []
    return {
      total: resultsData?.total || 0,
      qc: results.filter((r: Result) => r.result_type === 'qc_report').length,
      alignment: results.filter((r: Result) => r.result_type === 'alignment').length,
      quantification: results.filter((r: Result) => r.result_type === 'quantification').length,
      de: results.filter((r: Result) => r.result_type === 'de_analysis').length,
    }
  }, [resultsData])

  // Filter configurations
  const filterConfigs: FilterConfig[] = [
    {
      type: 'search',
      key: 'search',
      placeholder: 'Search results by type or ID...',
    },
    {
      type: 'select',
      key: 'result_type',
      label: 'Result Type',
      options: [
        { label: 'All Types', value: 'all' },
        { label: 'QC Report', value: 'qc_report' },
        { label: 'Alignment', value: 'alignment' },
        { label: 'Quantification', value: 'quantification' },
        { label: 'DE Analysis', value: 'de_analysis' },
      ],
    },
  ]

  // Result type tag colors and icons
  const getResultTypeConfig = (type: string) => {
    const configs: Record<string, { color: string; icon: React.ReactNode; label: string }> = {
      qc_report: {
        color: 'blue',
        icon: <CheckCircleOutlined />,
        label: 'QC Report',
      },
      alignment: {
        color: 'green',
        icon: <ThunderboltOutlined />,
        label: 'Alignment',
      },
      quantification: {
        color: 'purple',
        icon: <BarChartOutlined />,
        label: 'Quantification',
      },
      de_analysis: {
        color: 'orange',
        icon: <ExperimentOutlined />,
        label: 'DE Analysis',
      },
    }
    return configs[type] || { color: 'default', icon: null, label: type }
  }

  // Table columns
  const columns: ColumnsType<Result> = [
    {
      title: 'Result Type',
      dataIndex: 'result_type',
      key: 'result_type',
      width: 180,
      render: (type: string) => {
        const config = getResultTypeConfig(type)
        return (
          <Tag icon={config.icon} color={config.color}>
            {config.label}
          </Tag>
        )
      },
      filters: [
        { text: 'QC Report', value: 'qc_report' },
        { text: 'Alignment', value: 'alignment' },
        { text: 'Quantification', value: 'quantification' },
        { text: 'DE Analysis', value: 'de_analysis' },
      ],
      onFilter: (value, record) => record.result_type === value,
    },
    {
      title: 'Task ID',
      dataIndex: 'task_id',
      key: 'task_id',
      width: 280,
      ellipsis: true,
      render: (taskId: string) => (
        <Tooltip title={taskId}>
          <Text copyable={{ text: taskId }} style={{ fontSize: 12, fontFamily: 'monospace' }}>
            {taskId.slice(0, 8)}...{taskId.slice(-8)}
          </Text>
        </Tooltip>
      ),
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 180,
      sorter: (a, b) => dayjs(a.created_at).unix() - dayjs(b.created_at).unix(),
      render: (date: string) => (
        <Tooltip title={dayjs(date).format('YYYY-MM-DD HH:mm:ss')}>
          <Text type="secondary">{dayjs(date).fromNow()}</Text>
        </Tooltip>
      ),
    },
    {
      title: 'Status',
      dataIndex: 'metadata',
      key: 'status',
      width: 120,
      render: (metadata: any) => {
        const status = metadata?.status || 'unknown'
        const statusConfig: Record<string, { color: string; icon: React.ReactNode }> = {
          pass: { color: 'success', icon: <CheckCircleOutlined /> },
          warning: { color: 'warning', icon: <WarningOutlined /> },
          fail: { color: 'error', icon: <CloseCircleOutlined /> },
          unknown: { color: 'default', icon: null },
        }
        const config = statusConfig[status]
        return <StatusTag status={config.color} icon={config.icon} label={status} />
      },
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 150,
      fixed: 'right',
      render: (_, record) => {
        const items: MenuProps['items'] = [
          {
            key: 'view',
            label: 'View Details',
            icon: <EyeOutlined />,
            onClick: () => navigate(`/results/${record.id}`),
          },
          {
            key: 'visualize',
            label: 'Visualize',
            icon: <BarChartOutlined />,
            onClick: () => navigate(`/results/${record.id}#visualizations`),
          },
          {
            type: 'divider',
          },
          {
            key: 'download',
            label: 'Download',
            icon: <DownloadOutlined />,
            onClick: () => handleDownload(record),
          },
        ]

        return (
          <Space>
            <Button
              type="primary"
              size="small"
              icon={<EyeOutlined />}
              onClick={() => navigate(`/results/${record.id}`)}
            >
              View
            </Button>
            <Dropdown menu={{ items }} trigger={['click']}>
              <Button size="small" icon={<MoreOutlined />} />
            </Dropdown>
          </Space>
        )
      },
    },
  ]

  const handleDownload = async (result: Result) => {
    try {
      toast.info('Download functionality coming soon!')
      // TODO: Implement result download
    } catch (error: any) {
      toast.error(`Download failed: ${error.message}`)
    }
  }

  const rowSelection = {
    selectedRowKeys,
    onChange: (keys: React.Key[]) => setSelectedRowKeys(keys),
  }

  // Loading state
  if (loading && !resultsData) {
    return <PageSkeleton hasHeader rows={6} />
  }

  return (
    <div>
      {/* Header with Statistics */}
      <FadeIn direction="up" delay={0}>
        <Space direction="vertical" size="large" style={{ width: '100%', marginBottom: 24 }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
            <div>
              <Title level={2} style={{ margin: 0 }}>
                Analysis Results
              </Title>
              <Text type="secondary">Browse and visualize your analysis results</Text>
            </div>
            <Button
              type="primary"
              icon={<BarChartOutlined />}
              onClick={() => toast.info('Bulk visualization coming soon!')}
            >
              Visualize Selected
            </Button>
          </div>

          {/* Statistics Cards */}
          <StaggeredList staggerDelay={60} baseDelay={100} direction="up">
            <Row gutter={[16, 16]}>
              <Col xs={24} sm={12} md={6} lg={6}>
                <Card size="small">
                  <Statistic
                    title="Total Results"
                    value={stats.total}
                    prefix={<BarChartOutlined />}
                    valueStyle={{ color: '#1890ff' }}
                  />
                </Card>
              </Col>
              <Col xs={24} sm={12} md={6} lg={6}>
                <Card size="small">
                  <Statistic
                    title="QC Reports"
                    value={stats.qc}
                    prefix={<CheckCircleOutlined />}
                    valueStyle={{ color: '#52c41a' }}
                  />
                </Card>
              </Col>
              <Col xs={24} sm={12} md={6} lg={6}>
                <Card size="small">
                  <Statistic
                    title="Alignments"
                    value={stats.alignment}
                    prefix={<ThunderboltOutlined />}
                    valueStyle={{ color: '#722ed1' }}
                  />
                </Card>
              </Col>
              <Col xs={24} sm={12} md={6} lg={6}>
                <Card size="small">
                  <Statistic
                    title="DE Analysis"
                    value={stats.de}
                    prefix={<ExperimentOutlined />}
                    valueStyle={{ color: '#fa8c16' }}
                  />
                </Card>
              </Col>
            </Row>
          </StaggeredList>
        </Space>
      </FadeIn>

      {/* Filters */}
      <FadeIn direction="up" delay={100}>
        <FilterBar
          configs={filterConfigs}
          filters={filters}
          onFilterChange={(key, value) => {
            setFilter(key, value)
            resetPagination()
          }}
          onReset={() => {
            resetFilters()
            resetPagination()
          }}
          style={{ marginBottom: 16 }}
        />
      </FadeIn>

      {/* Results Table */}
      <FadeIn direction="up" delay={200}>
        <Card>
          <Table
            rowSelection={rowSelection}
            columns={columns}
            dataSource={filteredResults}
            rowKey="id"
            loading={loading}
            pagination={{
              current: pagination.current,
              pageSize: pagination.pageSize,
              total: resultsData?.total || 0,
              showSizeChanger: true,
              showQuickJumper: true,
              showTotal: (total, range) => `${range[0]}-${range[1]} of ${total} results`,
              onChange: handlePageChange,
            }}
            locale={{
              emptyText: (
                <EnhancedEmptyState
                  title="No Results Found"
                  description="No analysis results match your current filters. Try adjusting your search criteria or run a new analysis."
                  actionText="Reset Filters"
                  onAction={resetFilters}
                />
              ),
            }}
            scroll={{ x: 1000 }}
          />
        </Card>
      </FadeIn>

      {selectedRowKeys.length > 0 && (
        <div
          style={{
            position: 'fixed',
            bottom: 24,
            right: 24,
            background: '#fff',
            padding: '12px 24px',
            borderRadius: 8,
            boxShadow: '0 4px 12px rgba(0,0,0,0.15)',
            zIndex: 1000,
          }}
        >
          <Space>
            <Text strong>{selectedRowKeys.length} results selected</Text>
            <Button size="small" onClick={() => setSelectedRowKeys([])}>
              Clear
            </Button>
            <Button
              type="primary"
              size="small"
              icon={<DownloadOutlined />}
              onClick={() => toast.info('Bulk download coming soon!')}
            >
              Download Selected
            </Button>
          </Space>
        </div>
      )}
    </div>
  )
}

export default ResultList
