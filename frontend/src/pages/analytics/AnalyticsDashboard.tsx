import { useState, useEffect, useCallback } from 'react'
import { Card, Row, Col, Typography, Select, Space, Button, Statistic, Progress, Tag, Tabs, Alert } from 'antd'
import {
  LineChartOutlined,
  BarChartOutlined,
  PieChartOutlined,
  FundOutlined,
  RiseOutlined,
  FallOutlined,
  SwapOutlined,
  DownloadOutlined,
  ReloadOutlined,
  DashboardOutlined,
} from '@ant-design/icons'
import { PageHeader, StatisticCard, PageSkeleton, FadeIn, EnhancedEmptyState } from '@/components/common'
import type { AnalyticsSummary, TrendAnalysis } from '@/types/analytics'
import { analyticsService } from '@/services/analytics.service'
import { DesignTokens } from '@/styles/design-tokens'
import { logger } from '@/utils/logger'
import './AnalyticsDashboard.css'

const { Title, Text } = Typography
const { Option } = Select
const { TabPane } = Tabs

export const AnalyticsDashboard: React.FC = () => {
  const [loading, setLoading] = useState(false)
  const [initialLoad, setInitialLoad] = useState(true)
  const [period, setPeriod] = useState<'today' | 'week' | 'month' | 'year'>('week')
  const [summary, setSummary] = useState<AnalyticsSummary | null>(null)
  const [trends, setTrends] = useState<TrendAnalysis[]>([])

  const fetchAnalytics = useCallback(async () => {
    setLoading(true)
    try {
      const [summaryData, trendsData] = await Promise.all([
        analyticsService.getAnalyticsSummary(period),
        analyticsService.getTrendAnalysis(['Project Success Rate', 'Avg Quality Score', 'Pipeline Failures']),
      ])
      setSummary(summaryData)
      setTrends(trendsData)
    } catch (error) {
      logger.error('Failed to fetch analytics:', error)
    } finally {
      setLoading(false)
      setInitialLoad(false)
    }
  }, [period])

  useEffect(() => {
    fetchAnalytics()
  }, [fetchAnalytics])

  const getTrendIcon = (trend: 'up' | 'down' | 'stable') => {
    if (trend === 'up') {
      return <RiseOutlined style={{ color: DesignTokens.colors.success.main }} />
    }
    if (trend === 'down') {
      return <FallOutlined style={{ color: DesignTokens.colors.error.main }} />
    }
    return <SwapOutlined style={{ color: DesignTokens.colors.info.main }} />
  }

  const getTrendColor = (trend: 'up' | 'down' | 'stable', metric: string) => {
    // For some metrics, down is good (like failures)
    const reverseMetrics = ['failures', 'cost', 'pending']
    const isReverse = reverseMetrics.some((m) => metric.toLowerCase().includes(m))

    if (trend === 'up') {
      return isReverse ? DesignTokens.colors.error.main : DesignTokens.colors.success.main
    }
    if (trend === 'down') {
      return isReverse ? DesignTokens.colors.success.main : DesignTokens.colors.error.main
    }
    return DesignTokens.colors.info.main
  }

  // Show skeleton on initial load
  if (initialLoad && loading) {
    return <PageSkeleton hasHeader rows={8} />
  }

  if (!summary) {
    return null
  }

  return (
    <div className="analytics-dashboard">
      <FadeIn direction="up" delay={0} duration={300}>
        <PageHeader
          title="分析仪表板"
          subtitle="综合数据分析与可视化"
          icon={<DashboardOutlined />}
          extra={
            <Space>
              <Select value={period} onChange={setPeriod} style={{ width: 120 }}>
                <Option value="today">今天</Option>
                <Option value="week">本周</Option>
                <Option value="month">本月</Option>
                <Option value="year">本年</Option>
              </Select>
              <Button icon={<ReloadOutlined />} onClick={fetchAnalytics}>
                刷新
              </Button>
              <Button type="primary" icon={<DownloadOutlined />}>
                导出报告
              </Button>
            </Space>
          }
        />
      </FadeIn>

      {/* Trend Alerts */}
      {trends.length > 0 && (
        <FadeIn direction="up" delay={50} duration={300}>
          <Alert
            message="关键趋势"
            description={
              <Space direction="vertical" style={{ width: '100%' }} size="small">
                {trends.map((trend, idx) => (
                  <Text key={idx}>
                    {getTrendIcon(trend.trend)} <Text strong>{trend.metric}</Text>:{' '}
                    <Text style={{ color: getTrendColor(trend.trend, trend.metric) }}>
                      {trend.change > 0 ? '+' : ''}
                      {trend.change.toFixed(1)} ({trend.changePercent > 0 ? '+' : ''}
                      {trend.changePercent.toFixed(1)}%)
                    </Text>
                  </Text>
                ))}
              </Space>
            }
            type="info"
            showIcon
            closable
            style={{ marginBottom: 24 }}
          />
        </FadeIn>
      )}

      {/* Main Stats */}
      <FadeIn direction="up" delay={100} duration={300}>
        <Tabs defaultActiveKey="overview" size="large">
          <TabPane
            tab={
              <Space>
                <FundOutlined />
                总览
              </Space>
            }
            key="overview"
          >
            {/* Project Stats */}
            <Title level={4} style={{ marginTop: 0 }}>
              <LineChartOutlined /> 项目统计
            </Title>
            <Row gutter={[16, 16]} style={{ marginBottom: 32 }}>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="总项目数"
                  value={summary.projectStats.total}
                  icon={<FundOutlined />}
                  color={DesignTokens.colors.primary.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="进行中"
                  value={summary.projectStats.active}
                  icon={<LineChartOutlined />}
                  color={DesignTokens.colors.warning.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="已完成"
                  value={summary.projectStats.completed}
                  icon={<BarChartOutlined />}
                  color={DesignTokens.colors.success.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <Card className="stat-card">
                  <Statistic
                    title="成功率"
                    value={summary.projectStats.successRate}
                    suffix="%"
                    valueStyle={{ color: DesignTokens.colors.success.main }}
                  />
                  <Progress
                    percent={summary.projectStats.successRate}
                    strokeColor={DesignTokens.colors.success.main}
                    showInfo={false}
                    style={{ marginTop: 8 }}
                  />
                </Card>
              </Col>
            </Row>

            {/* Sample Stats */}
            <Title level={4}>
              <PieChartOutlined /> 样本统计
            </Title>
            <Row gutter={[16, 16]} style={{ marginBottom: 32 }}>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="总样本数"
                  value={summary.sampleStats.total}
                  icon={<PieChartOutlined />}
                  color={DesignTokens.colors.info.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="已处理"
                  value={summary.sampleStats.processed}
                  icon={<BarChartOutlined />}
                  color={DesignTokens.colors.success.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <StatisticCard
                  title="待处理"
                  value={summary.sampleStats.pending}
                  icon={<LineChartOutlined />}
                  color={DesignTokens.colors.warning.main}
                />
              </Col>
              <Col xs={12} sm={6}>
                <Card className="stat-card">
                  <Statistic
                    title="平均质量分"
                    value={summary.sampleStats.avgQualityScore}
                    suffix="/ 100"
                    valueStyle={{
                      color:
                        summary.sampleStats.avgQualityScore >= 80
                          ? DesignTokens.colors.success.main
                          : DesignTokens.colors.warning.main,
                    }}
                  />
                  <Progress
                    percent={summary.sampleStats.avgQualityScore}
                    strokeColor={
                      summary.sampleStats.avgQualityScore >= 80
                        ? DesignTokens.colors.success.main
                        : DesignTokens.colors.warning.main
                    }
                    showInfo={false}
                    style={{ marginTop: 8 }}
                  />
                </Card>
              </Col>
            </Row>

            {/* Pipeline Stats */}
            <Title level={4}>
              <FundOutlined /> 流程统计
            </Title>
            <Row gutter={[16, 16]}>
              <Col xs={24} md={12} lg={6}>
                <Card className="stat-card">
                  <Statistic title="总运行次数" value={summary.pipelineStats.totalRuns} prefix={<FundOutlined />} />
                </Card>
              </Col>
              <Col xs={24} md={12} lg={6}>
                <Card className="stat-card">
                  <Statistic
                    title="成功运行"
                    value={summary.pipelineStats.successfulRuns}
                    valueStyle={{ color: DesignTokens.colors.success.main }}
                  />
                </Card>
              </Col>
              <Col xs={24} md={12} lg={6}>
                <Card className="stat-card">
                  <Statistic
                    title="平均运行时间"
                    value={summary.pipelineStats.avgRuntime}
                    valueStyle={{ fontSize: 20 }}
                  />
                </Card>
              </Col>
              <Col xs={24} md={12} lg={6}>
                <Card className="stat-card">
                  <Statistic
                    title="平均成本"
                    value={summary.pipelineStats.avgCost}
                    prefix="$"
                    precision={2}
                    valueStyle={{ fontSize: 20 }}
                  />
                </Card>
              </Col>
            </Row>
          </TabPane>

          <TabPane
            tab={
              <Space>
                <BarChartOutlined />
                存储分析
              </Space>
            }
            key="storage"
          >
            <Row gutter={[24, 24]}>
              <Col xs={24} lg={12}>
                <Card title="存储使用情况" className="storage-card">
                  <Space direction="vertical" style={{ width: '100%' }} size="large">
                    <div>
                      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 8 }}>
                        <Text>总使用量</Text>
                        <Text strong>
                          {summary.storageStats.totalUsed} GB / {summary.storageStats.totalAvailable} GB
                        </Text>
                      </div>
                      <Progress
                        percent={(summary.storageStats.totalUsed / summary.storageStats.totalAvailable) * 100}
                        strokeColor={
                          summary.storageStats.totalUsed / summary.storageStats.totalAvailable > 0.8
                            ? DesignTokens.colors.error.main
                            : DesignTokens.colors.primary.main
                        }
                      />
                    </div>

                    <div>
                      <Title level={5}>按类型分布</Title>
                      {Object.entries(summary.storageStats.usageByType).map(([type, size]) => (
                        <div key={type} style={{ marginBottom: 12 }}>
                          <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
                            <Text>{type.toUpperCase()}</Text>
                            <Text>{size} GB</Text>
                          </div>
                          <Progress percent={(size / summary.storageStats.totalUsed) * 100} showInfo={false} />
                        </div>
                      ))}
                    </div>

                    <div>
                      <Tag
                        color={
                          summary.storageStats.trend === 'increasing'
                            ? 'orange'
                            : summary.storageStats.trend === 'decreasing'
                              ? 'green'
                              : 'blue'
                        }
                      >
                        趋势:{' '}
                        {summary.storageStats.trend === 'increasing'
                          ? '增长中'
                          : summary.storageStats.trend === 'decreasing'
                            ? '下降中'
                            : '稳定'}
                      </Tag>
                    </div>
                  </Space>
                </Card>
              </Col>

              <Col xs={24} lg={12}>
                <Card title="活动统计" className="activity-card">
                  <Space direction="vertical" style={{ width: '100%' }} size="large">
                    <Statistic title="日活跃用户" value={summary.activityStats.dailyActiveUsers} suffix="人" />
                    <Statistic title="总会话数" value={summary.activityStats.totalSessions} />
                    <Statistic title="平均会话时长" value={summary.activityStats.avgSessionDuration} />
                    <div>
                      <Text type="secondary">使用高峰时段: </Text>
                      <Text strong>{summary.activityStats.peakUsageHour}:00</Text>
                    </div>
                  </Space>
                </Card>
              </Col>
            </Row>
          </TabPane>

          <TabPane
            tab={
              <Space>
                <LineChartOutlined />
                趋势分析
              </Space>
            }
            key="trends"
          >
            <EnhancedEmptyState
              type="noData"
              title="Trend charts coming soon"
              description="趋势图表将在此显示"
              action={{
                text: '查看详细趋势',
                onClick: () => {},
              }}
              size="default"
            />
          </TabPane>
        </Tabs>
      </FadeIn>
    </div>
  )
}
