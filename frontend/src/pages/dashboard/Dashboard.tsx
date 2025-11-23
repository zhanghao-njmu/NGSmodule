/**
 * Dashboard Page
 * Refactored: Unified color usage with CSS variables
 */
import { Card, Row, Col, Statistic, Typography, Space, Tag, Progress } from 'antd'
import {
  ProjectOutlined,
  ExperimentOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  DatabaseOutlined,
} from '@ant-design/icons'
import { authStore } from '@/store/authStore'
import { PageSkeleton, FadeIn, StaggeredList, EnhancedEmptyState } from '@/components/common'
import { useAsync } from '@/hooks'
import { statsService } from '@/services/stats.service'
import styles from './Dashboard.module.css'

const { Title, Text } = Typography

export const Dashboard: React.FC = () => {
  const { user } = authStore()

  // Using useAsync hook eliminates 20+ lines of boilerplate
  const {
    data: stats,
    loading,
    error,
    execute: loadStats,
  } = useAsync(
    async () => {
      // 获取仪表板统计数据并添加用户存储信息
      const dashboardStats = await statsService.getDashboardStats()
      return {
        ...dashboardStats,
        storageUsed: user?.storage_used || 0,
        storageQuota: user?.storage_quota || 107374182400, // 100GB 默认配额
      }
    },
    {
      immediate: true,
      onError: (err) => console.error('Failed to load dashboard stats:', err),
    },
  )

  const storagePercent = stats ? Math.round((stats.storageUsed / stats.storageQuota) * 100) : 0

  const formatBytes = (bytes: number) => {
    if (bytes === 0) {
      return '0 Bytes'
    }
    const k = 1024
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return Math.round((bytes / Math.pow(k, i)) * 100) / 100 + ' ' + sizes[i]
  }

  // Loading 状态 - Use PageSkeleton
  if (loading) {
    return <PageSkeleton hasHeader rows={4} />
  }

  // 错误状态
  if (error) {
    return (
      <div className={styles.dashboard}>
        <FadeIn>
          <EnhancedEmptyState
            type="error"
            title="Error Loading Dashboard"
            description={error?.message || 'Failed to load statistics. Please try again later.'}
            action={{
              text: 'Retry',
              onClick: loadStats,
            }}
            size="default"
          />
        </FadeIn>
      </div>
    )
  }

  // 无数据状态
  if (!stats) {
    return (
      <div className={styles.dashboard}>
        <FadeIn>
          <EnhancedEmptyState
            type="noData"
            title="No data available"
            description="Dashboard statistics are not available at this moment"
            size="default"
          />
        </FadeIn>
      </div>
    )
  }

  return (
    <div className={styles.dashboard}>
      {/* Header with Fade In */}
      <FadeIn direction="up" delay={0} duration={300}>
        <div className={styles.header}>
          <Space direction="vertical" size={4}>
            <Title level={2} style={{ margin: 0 }}>
              Welcome back, {user?.full_name || user?.username}! 👋
            </Title>
            <Text type="secondary">Here&apos;s an overview of your bioinformatics workspace</Text>
          </Space>
          <Tag color={user?.role === 'admin' ? 'gold' : 'blue'}>
            {user?.role === 'admin' ? 'Administrator' : 'User'}
          </Tag>
        </div>
      </FadeIn>

      {/* Statistics Cards with Staggered Animation */}
      <StaggeredList staggerDelay={80} baseDelay={100} direction="up">
        <Row gutter={[24, 24]}>
          <Col xs={24} sm={12} lg={6}>
            <Card bordered={false} className={styles.statCard}>
              <Statistic
                title="Total Projects"
                value={stats.totalProjects}
                prefix={<ProjectOutlined />}
                valueStyle={{ color: 'var(--color-primary)' }}
              />
            </Card>
          </Col>

          <Col xs={24} sm={12} lg={6}>
            <Card bordered={false} className={styles.statCard}>
              <Statistic
                title="Running Tasks"
                value={stats.runningTasks}
                prefix={<ClockCircleOutlined />}
                valueStyle={{ color: 'var(--color-warning)' }}
              />
            </Card>
          </Col>

          <Col xs={24} sm={12} lg={6}>
            <Card bordered={false} className={styles.statCard}>
              <Statistic
                title="Completed Tasks"
                value={stats.completedTasks}
                prefix={<CheckCircleOutlined />}
                valueStyle={{ color: 'var(--color-success)' }}
              />
            </Card>
          </Col>

          <Col xs={24} sm={12} lg={6}>
            <Card bordered={false} className={styles.statCard}>
              <Statistic
                title="Storage Used"
                value={storagePercent}
                suffix="%"
                prefix={<DatabaseOutlined />}
                valueStyle={{ color: storagePercent > 80 ? 'var(--color-error)' : 'var(--color-primary)' }}
              />
              <div style={{ marginTop: 12 }}>
                <Progress
                  percent={storagePercent}
                  strokeColor={storagePercent > 80 ? 'var(--color-error)' : 'var(--color-primary)'}
                  showInfo={false}
                  size="small"
                />
                <Text type="secondary" style={{ fontSize: 12, marginTop: 4, display: 'block' }}>
                  {formatBytes(stats.storageUsed)} / {formatBytes(stats.storageQuota)}
                </Text>
              </div>
            </Card>
          </Col>
        </Row>
      </StaggeredList>

      {/* Content Area with Fade In */}
      <FadeIn direction="up" delay={400} duration={300}>
        <Row gutter={[24, 24]} style={{ marginTop: 24 }}>
          <Col xs={24} lg={16}>
            <Card
              title="Recent Projects"
              bordered={false}
              extra={<a href="/items">View All</a>}
              className={styles.card}
            >
              <EnhancedEmptyState
                type="noData"
                title="No items yet"
                description="Create your first project to start analyzing NGS data"
                action={{
                  text: 'Create Project',
                  onClick: () => (window.location.href = '/items'),
                }}
                size="small"
              />
            </Card>
          </Col>

          <Col xs={24} lg={8}>
            <Card title="Quick Actions" bordered={false} className={styles.card}>
              <Space direction="vertical" style={{ width: '100%' }} size="middle">
                <Card.Grid
                  hoverable
                  style={{ width: '100%', cursor: 'pointer' }}
                  onClick={() => (window.location.href = '/items')}
                >
                  <Space>
                    <ProjectOutlined style={{ fontSize: 20, color: 'var(--color-primary)' }} />
                    <div>
                      <Text strong>Create Project</Text>
                      <br />
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        Start a new analysis project
                      </Text>
                    </div>
                  </Space>
                </Card.Grid>

                <Card.Grid hoverable style={{ width: '100%', cursor: 'pointer' }}>
                  <Space>
                    <ExperimentOutlined style={{ fontSize: 20, color: 'var(--color-success)' }} />
                    <div>
                      <Text strong>Run Pipeline</Text>
                      <br />
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        Execute NGS analysis workflow
                      </Text>
                    </div>
                  </Space>
                </Card.Grid>

                <Card.Grid hoverable style={{ width: '100%', cursor: 'pointer' }}>
                  <Space>
                    <DatabaseOutlined style={{ fontSize: 20, color: 'var(--color-warning)' }} />
                    <div>
                      <Text strong>Upload Data</Text>
                      <br />
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        Upload FASTQ or BAM files
                      </Text>
                    </div>
                  </Space>
                </Card.Grid>
              </Space>
            </Card>
          </Col>
        </Row>
      </FadeIn>
    </div>
  )
}

export default Dashboard
