import { useState, useEffect } from 'react'
import { Card, Row, Col, Typography, Statistic, Progress, Tag, Space, Button, Alert, List, Badge, Spin } from 'antd'
import {
  CheckCircleOutlined,
  WarningOutlined,
  CloseCircleOutlined,
  ReloadOutlined,
  DesktopOutlined,
  DatabaseOutlined,
  CloudOutlined,
  ApiOutlined,
} from '@ant-design/icons'
import { DesignTokens } from '@/styles/design-tokens'
import type { SystemHealth, SystemMetrics, SystemAlert } from '@/types/admin'
import './SystemMonitor.css'

const { Text } = Typography

export const SystemMonitor: React.FC = () => {
  const [loading, setLoading] = useState(false)
  const [health, setHealth] = useState<SystemHealth | null>(null)
  const [metrics, setMetrics] = useState<SystemMetrics | null>(null)
  const [alerts, setAlerts] = useState<SystemAlert[]>([])

  useEffect(() => {
    fetchSystemData()
    const interval = setInterval(fetchSystemData, 30000) // Refresh every 30s
    return () => clearInterval(interval)
  }, [])

  const fetchSystemData = async () => {
    setLoading(true)
    try {
      // TODO: Replace with actual API calls
      // const healthData = await enhancedAdminService.getSystemHealth()
      // const metricsData = await enhancedAdminService.getSystemMetrics()
      // const alertsData = await enhancedAdminService.getAlerts()

      // Mock data
      const mockHealth: SystemHealth = {
        overall: 'healthy',
        components: {
          api: { status: 'operational', latency: 45 },
          database: { status: 'operational', latency: 12 },
          storage: { status: 'operational' },
          queue: { status: 'operational' },
        },
        metrics: {
          cpu: 45.2,
          memory: 62.8,
          disk: 58.3,
        },
      }

      const mockMetrics: SystemMetrics = {
        cpu: { usage: 45.2, load: [1.2, 1.5, 1.8] },
        memory: { used: 6280, total: 10000, usagePercent: 62.8 },
        disk: { used: 583, total: 1000, usagePercent: 58.3 },
      }

      const mockAlerts: SystemAlert[] = [
        {
          id: '1',
          type: 'warning',
          severity: 'medium',
          title: 'High Memory Usage',
          message: 'Memory usage above 60%',
          timestamp: new Date().toISOString(),
          resolved: false,
        },
      ]

      setHealth(mockHealth)
      setMetrics(mockMetrics)
      setAlerts(mockAlerts)
    } catch (error) {
      console.error('Failed to fetch system data:', error)
    } finally {
      setLoading(false)
    }
  }

  const getHealthColor = (status: string): string => {
    if (status === 'healthy' || status === 'operational') {
      return DesignTokens.colors.success.main
    }
    if (status === 'degraded') {
      return DesignTokens.colors.warning.main
    }
    return DesignTokens.colors.error.main
  }

  const getHealthIcon = (status: string) => {
    if (status === 'healthy' || status === 'operational') {
      return <CheckCircleOutlined />
    }
    if (status === 'degraded') {
      return <WarningOutlined />
    }
    return <CloseCircleOutlined />
  }

  const getSeverityColor = (severity: string): string => {
    const colors: Record<string, string> = {
      critical: 'red',
      high: 'orange',
      medium: 'gold',
      low: 'blue',
    }
    return colors[severity] || 'default'
  }

  if (loading && !health) {
    return (
      <div style={{ textAlign: 'center', padding: '48px' }}>
        <Spin size="large" />
        <Text type="secondary" style={{ display: 'block', marginTop: 16 }}>
          Loading system status...
        </Text>
      </div>
    )
  }

  if (!health || !metrics) {
    return null
  }

  return (
    <div className="system-monitor">
      {/* Overall Status */}
      <Card
        title={
          <Space>
            {getHealthIcon(health.overall)}
            <Text strong>System Health</Text>
          </Space>
        }
        extra={
          <Space>
            <Tag color={getHealthColor(health.overall)}>{health.overall.toUpperCase()}</Tag>
            <Button icon={<ReloadOutlined />} onClick={fetchSystemData} loading={loading}>
              Refresh
            </Button>
          </Space>
        }
        className="health-card"
      >
        <Row gutter={[16, 16]}>
          {/* Component Status */}
          {Object.entries(health.components).map(([name, component]) => (
            <Col xs={12} md={6} key={name}>
              <Card className="component-card" size="small">
                <Space direction="vertical" style={{ width: '100%' }} align="center">
                  <div style={{ fontSize: 24, color: getHealthColor(component.status) }}>
                    {name === 'api' && <ApiOutlined />}
                    {name === 'database' && <DatabaseOutlined />}
                    {name === 'storage' && <CloudOutlined />}
                    {name === 'queue' && <DesktopOutlined />}
                  </div>
                  <Text strong style={{ textTransform: 'capitalize' }}>
                    {name}
                  </Text>
                  <Tag color={getHealthColor(component.status)}>{component.status}</Tag>
                  {component.latency && (
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      {component.latency}ms
                    </Text>
                  )}
                </Space>
              </Card>
            </Col>
          ))}
        </Row>
      </Card>

      {/* Resource Metrics */}
      <Row gutter={[16, 16]} style={{ marginTop: 16 }}>
        <Col xs={24} md={8}>
          <Card className="metric-card">
            <Statistic
              title="CPU Usage"
              value={metrics.cpu.usage}
              suffix="%"
              valueStyle={{
                color: metrics.cpu.usage > 80 ? DesignTokens.colors.error.main : DesignTokens.colors.success.main,
              }}
            />
            <Progress
              percent={metrics.cpu.usage}
              strokeColor={
                metrics.cpu.usage > 80
                  ? DesignTokens.colors.error.main
                  : metrics.cpu.usage > 60
                    ? DesignTokens.colors.warning.main
                    : DesignTokens.colors.success.main
              }
              style={{ marginTop: 12 }}
            />
            <Text type="secondary" style={{ fontSize: 12 }}>
              Load: {metrics.cpu.load.join(', ')}
            </Text>
          </Card>
        </Col>

        <Col xs={24} md={8}>
          <Card className="metric-card">
            <Statistic
              title="Memory Usage"
              value={metrics.memory.usagePercent}
              suffix="%"
              valueStyle={{
                color:
                  metrics.memory.usagePercent > 80 ? DesignTokens.colors.error.main : DesignTokens.colors.success.main,
              }}
            />
            <Progress
              percent={metrics.memory.usagePercent}
              strokeColor={
                metrics.memory.usagePercent > 80
                  ? DesignTokens.colors.error.main
                  : metrics.memory.usagePercent > 60
                    ? DesignTokens.colors.warning.main
                    : DesignTokens.colors.success.main
              }
              style={{ marginTop: 12 }}
            />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {metrics.memory.used} MB / {metrics.memory.total} MB
            </Text>
          </Card>
        </Col>

        <Col xs={24} md={8}>
          <Card className="metric-card">
            <Statistic
              title="Disk Usage"
              value={metrics.disk.usagePercent}
              suffix="%"
              valueStyle={{
                color:
                  metrics.disk.usagePercent > 80 ? DesignTokens.colors.error.main : DesignTokens.colors.success.main,
              }}
            />
            <Progress
              percent={metrics.disk.usagePercent}
              strokeColor={
                metrics.disk.usagePercent > 80
                  ? DesignTokens.colors.error.main
                  : metrics.disk.usagePercent > 60
                    ? DesignTokens.colors.warning.main
                    : DesignTokens.colors.success.main
              }
              style={{ marginTop: 12 }}
            />
            <Text type="secondary" style={{ fontSize: 12 }}>
              {metrics.disk.used} GB / {metrics.disk.total} GB
            </Text>
          </Card>
        </Col>
      </Row>

      {/* Active Alerts */}
      {alerts.length > 0 && (
        <Card
          title={
            <Space>
              <Badge count={alerts.length} />
              <Text strong>Active Alerts</Text>
            </Space>
          }
          style={{ marginTop: 16 }}
        >
          <List
            dataSource={alerts}
            renderItem={(alert) => (
              <Alert
                type={alert.type}
                message={
                  <Space style={{ width: '100%', justifyContent: 'space-between' }}>
                    <Space>
                      <Tag color={getSeverityColor(alert.severity)}>{alert.severity.toUpperCase()}</Tag>
                      <Text strong>{alert.title}</Text>
                    </Space>
                    <Text type="secondary" style={{ fontSize: 12 }}>
                      {new Date(alert.timestamp).toLocaleString()}
                    </Text>
                  </Space>
                }
                description={alert.message}
                showIcon
                closable
                onClose={() => {
                  // TODO: Resolve alert
                  setAlerts(alerts.filter((a) => a.id !== alert.id))
                }}
                style={{ marginBottom: 8 }}
              />
            )}
          />
        </Card>
      )}
    </div>
  )
}
