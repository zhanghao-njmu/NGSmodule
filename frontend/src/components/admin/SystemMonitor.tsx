import { useState, useEffect, useCallback } from 'react'
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
import { enhancedAdminService } from '@/services/admin.enhanced.service'
import type { SystemHealth, SystemMetrics, SystemAlert } from '@/types/admin'
import { logger } from '@/utils/logger'
import './SystemMonitor.css'

const { Text } = Typography

export const SystemMonitor: React.FC = () => {
  const [loading, setLoading] = useState(false)
  const [health, setHealth] = useState<SystemHealth | null>(null)
  const [metrics, setMetrics] = useState<SystemMetrics | null>(null)
  const [alerts, setAlerts] = useState<SystemAlert[]>([])

  const fetchSystemData = useCallback(async () => {
    setLoading(true)
    try {
      const [healthData, metricsData, alertsData] = await Promise.all([
        enhancedAdminService.getSystemHealth(),
        enhancedAdminService.getSystemMetrics(),
        enhancedAdminService.getAlerts(false),
      ])

      setHealth(healthData as SystemHealth)
      setMetrics(metricsData as SystemMetrics)
      setAlerts(alertsData as SystemAlert[])
    } catch (error) {
      logger.error('Failed to fetch system data:', error)
    } finally {
      setLoading(false)
    }
  }, [])

  useEffect(() => {
    fetchSystemData()
    const interval = setInterval(fetchSystemData, 30000) // Refresh every 30s
    return () => clearInterval(interval)
  }, [fetchSystemData])

  const handleResolveAlert = async (alertId: string) => {
    try {
      await enhancedAdminService.resolveAlert(alertId)
      setAlerts(alerts.filter((a) => a.id !== alertId))
    } catch (error) {
      logger.error('Failed to resolve alert:', error)
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
                onClose={() => handleResolveAlert(alert.id)}
                style={{ marginBottom: 8 }}
              />
            )}
          />
        </Card>
      )}
    </div>
  )
}
