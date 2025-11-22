/**
 * Dashboard Page
 */
import React, { useEffect, useState } from 'react'
import { Card, Row, Col, Statistic, Typography, Space, Tag, Spin, Alert } from 'antd'
import {
  ProjectOutlined,
  ExperimentOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
  DatabaseOutlined,
} from '@ant-design/icons'
import { authStore } from '@/store/authStore'
import statsService, { type DashboardStats } from '@/services/stats.service'
import styles from './Dashboard.module.css'

const { Title, Text, Paragraph } = Typography

export const Dashboard: React.FC = () => {
  const { user } = authStore()
  const [stats, setStats] = useState<DashboardStats | null>(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  // 加载统计数据
  useEffect(() => {
    loadStats()
  }, [])

  const loadStats = async () => {
    try {
      setLoading(true)
      setError(null)
      const data = await statsService.getDashboardStats()

      // 将用户存储信息添加到统计数据
      setStats({
        ...data,
        storageUsed: user?.storage_used || 0,
        storageQuota: user?.storage_quota || 107374182400, // 100GB 默认配额
      })
    } catch (err) {
      console.error('Failed to load dashboard stats:', err)
      setError('Failed to load statistics. Please try again later.')
    } finally {
      setLoading(false)
    }
  }

  const storagePercent = stats
    ? Math.round((stats.storageUsed / stats.storageQuota) * 100)
    : 0

  const formatBytes = (bytes: number) => {
    if (bytes === 0) return '0 Bytes'
    const k = 1024
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return Math.round((bytes / Math.pow(k, i)) * 100) / 100 + ' ' + sizes[i]
  }

  // Loading 状态
  if (loading) {
    return (
      <div className={styles.dashboard}>
        <div style={{ textAlign: 'center', padding: '100px 0' }}>
          <Spin size="large" tip="Loading dashboard..." />
        </div>
      </div>
    )
  }

  // 错误状态
  if (error) {
    return (
      <div className={styles.dashboard}>
        <Alert
          message="Error Loading Dashboard"
          description={error}
          type="error"
          showIcon
          action={
            <a onClick={loadStats} style={{ cursor: 'pointer' }}>
              Retry
            </a>
          }
        />
      </div>
    )
  }

  // 无数据状态
  if (!stats) {
    return (
      <div className={styles.dashboard}>
        <Alert message="No data available" type="info" />
      </div>
    )
  }

  return (
    <div className={styles.dashboard}>
      <div className={styles.header}>
        <Space direction="vertical" size={4}>
          <Title level={2} style={{ margin: 0 }}>
            Welcome back, {user?.full_name || user?.username}! 👋
          </Title>
          <Text type="secondary">
            Here's an overview of your bioinformatics workspace
          </Text>
        </Space>
        <Tag color={user?.role === 'admin' ? 'gold' : 'blue'}>
          {user?.role === 'admin' ? 'Administrator' : 'User'}
        </Tag>
      </div>

      <Row gutter={[24, 24]}>
        <Col xs={24} sm={12} lg={6}>
          <Card bordered={false} className={styles.statCard}>
            <Statistic
              title="Total Projects"
              value={stats.totalProjects}
              prefix={<ProjectOutlined />}
              valueStyle={{ color: '#2196F3' }}
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card bordered={false} className={styles.statCard}>
            <Statistic
              title="Running Tasks"
              value={stats.runningTasks}
              prefix={<ClockCircleOutlined />}
              valueStyle={{ color: '#FF9800' }}
            />
          </Card>
        </Col>

        <Col xs={24} sm={12} lg={6}>
          <Card bordered={false} className={styles.statCard}>
            <Statistic
              title="Completed Tasks"
              value={stats.completedTasks}
              prefix={<CheckCircleOutlined />}
              valueStyle={{ color: '#4CAF50' }}
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
              valueStyle={{ color: storagePercent > 80 ? '#F44336' : '#2196F3' }}
            />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {formatBytes(stats.storageUsed)} / {formatBytes(stats.storageQuota)}
            </Text>
          </Card>
        </Col>
      </Row>

      <Row gutter={[24, 24]} style={{ marginTop: 24 }}>
        <Col xs={24} lg={16}>
          <Card
            title="Recent Projects"
            bordered={false}
            extra={<a href="/projects">View All</a>}
            className={styles.card}
          >
            <div className={styles.emptyState}>
              <ExperimentOutlined style={{ fontSize: 48, color: '#ccc' }} />
              <Title level={4} type="secondary">
                No projects yet
              </Title>
              <Paragraph type="secondary">
                Create your first project to start analyzing NGS data
              </Paragraph>
            </div>
          </Card>
        </Col>

        <Col xs={24} lg={8}>
          <Card
            title="Quick Actions"
            bordered={false}
            className={styles.card}
          >
            <Space direction="vertical" style={{ width: '100%' }} size="middle">
              <Card.Grid
                hoverable
                style={{ width: '100%', cursor: 'pointer' }}
                onClick={() => window.location.href = '/projects'}
              >
                <Space>
                  <ProjectOutlined style={{ fontSize: 20, color: '#2196F3' }} />
                  <div>
                    <Text strong>Create Project</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      Start a new analysis project
                    </Text>
                  </div>
                </Space>
              </Card.Grid>

              <Card.Grid
                hoverable
                style={{ width: '100%', cursor: 'pointer' }}
              >
                <Space>
                  <ExperimentOutlined style={{ fontSize: 20, color: '#4CAF50' }} />
                  <div>
                    <Text strong>Run Pipeline</Text>
                    <br />
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      Execute NGS analysis workflow
                    </Text>
                  </div>
                </Space>
              </Card.Grid>

              <Card.Grid
                hoverable
                style={{ width: '100%', cursor: 'pointer' }}
              >
                <Space>
                  <DatabaseOutlined style={{ fontSize: 20, color: '#FF9800' }} />
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
    </div>
  )
}

export default Dashboard
