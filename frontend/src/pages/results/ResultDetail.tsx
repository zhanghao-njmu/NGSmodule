/**
 * Result Detail Page - Analysis results visualization
 */
import type React from 'react'
import { useEffect } from 'react'
import { useParams, useNavigate } from 'react-router-dom'
import { Card, Row, Col, Statistic, Tabs, Alert, Button, Space, Tag, Descriptions } from 'antd'
import { ArrowLeftOutlined, CheckCircleOutlined, WarningOutlined, CloseCircleOutlined } from '@ant-design/icons'
import resultService from '@/services/result.service'
import { PageSkeleton, FadeIn } from '@/components/common'
import { LineChart, BarChart, PieChart, ScatterPlot } from '@/components/charts'
import { useAsync } from '@/hooks'
import type { ResultVisualizationData, ChartData } from '@/types/result'

export const ResultDetail: React.FC = () => {
  const { id } = useParams<{ id: string }>()
  const navigate = useNavigate()

  // Using useAsync hook eliminates 15+ lines of boilerplate
  const {
    data: vizData,
    loading,
    error,
    execute: loadVisualizationData,
  } = useAsync(() => resultService.getVisualizationData(id!), { immediate: false })

  useEffect(() => {
    if (id) {
      loadVisualizationData()
    }
  }, [id])

  if (loading) {
    return <PageSkeleton hasHeader rows={6} />
  }

  if (error || !vizData) {
    return (
      <FadeIn>
        <Alert
          message="Error Loading Results"
          description={error?.message || 'No data available'}
          type="error"
          showIcon
          action={<Button onClick={loadVisualizationData}>Retry</Button>}
        />
      </FadeIn>
    )
  }

  // Status icon and color
  const statusConfig = {
    pass: { icon: <CheckCircleOutlined />, color: 'success' as const, text: 'Pass' },
    warning: { icon: <WarningOutlined />, color: 'warning' as const, text: 'Warning' },
    fail: { icon: <CloseCircleOutlined />, color: 'error' as const, text: 'Fail' },
  }

  const status = statusConfig[vizData.status || 'pass']

  return (
    <div>
      {/* Header */}
      <FadeIn direction="up" delay={0}>
        <Card bordered={false} style={{ marginBottom: 24 }}>
          <Space direction="vertical" size="large" style={{ width: '100%' }}>
            <Space>
              <Button icon={<ArrowLeftOutlined />} onClick={() => navigate(-1)}>
                Back
              </Button>
              <h2 style={{ margin: 0 }}>{vizData.type.replace('_', ' ').toUpperCase()} Results</h2>
              {vizData.status && (
                <Tag icon={status.icon} color={status.color}>
                  {status.text}
                </Tag>
              )}
            </Space>

            {/* Key Metrics */}
            <Row gutter={[16, 16]}>
              {Object.entries(vizData.metrics).map(([key, value]) => (
                <Col xs={24} sm={12} md={6} key={key}>
                  <Card size="small">
                    <Statistic
                      title={key.replace(/_/g, ' ').toUpperCase()}
                      value={typeof value === 'number' ? value : String(value)}
                      precision={typeof value === 'number' && value % 1 !== 0 ? 2 : 0}
                    />
                  </Card>
                </Col>
              ))}
            </Row>
          </Space>
        </Card>
      </FadeIn>

      {/* Charts */}
      <FadeIn direction="up" delay={100}>
        <Tabs
          defaultActiveKey="charts"
          items={[
            {
              key: 'charts',
              label: 'Visualizations',
              children: (
                <Row gutter={[16, 16]}>
                  {Object.entries(vizData.charts).map(([chartName, chartData]) => (
                    <Col xs={24} lg={12} key={chartName}>
                      {renderChart(chartName, chartData as ChartData)}
                    </Col>
                  ))}
                </Row>
              ),
            },
            ...(vizData.top_genes || vizData.significant_genes
              ? [
                  {
                    key: 'genes',
                    label: 'Gene List',
                    children: renderGeneTable(vizData),
                  },
                ]
              : []),
            {
              key: 'raw',
              label: 'Raw Data',
              children: (
                <Card>
                  <pre style={{ maxHeight: 600, overflow: 'auto' }}>{JSON.stringify(vizData, null, 2)}</pre>
                </Card>
              ),
            },
          ]}
        />
      </FadeIn>
    </div>
  )
}

function renderChart(name: string, data: ChartData) {
  const title = name.replace(/_/g, ' ').toUpperCase()

  switch (data.type) {
    case 'line':
      return <LineChart data={[{ x: data.x || [], y: data.y || [], name: title }]} title={title} smooth />

    case 'bar':
      return (
        <BarChart
          data={[{ categories: data.categories || data.x || [], values: data.values || data.y || [] }]}
          title={title}
        />
      )

    case 'pie':
      return (
        <PieChart
          data={
            data.categories && data.values
              ? data.categories.map((cat, i) => ({ name: cat, value: data.values![i] }))
              : []
          }
          title={title}
        />
      )

    case 'scatter':
      return (
        <ScatterPlot
          data={[
            {
              x: data.x || [],
              y: data.y || [],
              labels: data.labels,
              name: title,
            },
          ]}
          title={title}
        />
      )

    case 'area':
      return <LineChart data={[{ x: data.x || [], y: data.y || [], name: title }]} title={title} smooth showArea />

    case 'histogram':
      // For histogram, we can use bar chart
      return <BarChart data={[{ categories: data.x?.map(String) || [], values: data.y || [] }]} title={title} />

    default:
      return <Card title={title}>Unsupported chart type: {data.type}</Card>
  }
}

function renderGeneTable(vizData: ResultVisualizationData) {
  const genes = vizData.significant_genes || vizData.top_genes || []

  return (
    <Card>
      <Descriptions title="Gene Information" column={1} bordered size="small">
        {genes.slice(0, 20).map((gene, index) => (
          <Descriptions.Item key={index} label={gene.gene}>
            {gene.expression !== undefined && `Expression: ${gene.expression.toFixed(2)} | `}
            {gene.log2_fold_change !== undefined && `Log2FC: ${gene.log2_fold_change.toFixed(2)} | `}
            {gene.p_value !== undefined && `P-value: ${gene.p_value.toExponential(2)}`}
            {gene.significant && (
              <Tag color="red" style={{ marginLeft: 8 }}>
                Significant
              </Tag>
            )}
          </Descriptions.Item>
        ))}
      </Descriptions>
      {genes.length > 20 && (
        <div style={{ marginTop: 16, textAlign: 'center' }}>
          <Tag>Showing top 20 of {genes.length} genes</Tag>
        </div>
      )}
    </Card>
  )
}

export default ResultDetail
