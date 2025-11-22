import { useState, useEffect } from 'react'
import {
  Card,
  Space,
  Typography,
  Tag,
  Button,
  Progress,
  Alert,
  Spin,
  Collapse,
  Statistic,
  Row,
  Col,
  Tooltip,
} from 'antd'
import {
  CheckCircleOutlined,
  WarningOutlined,
  CloseCircleOutlined,
  InfoCircleOutlined,
  ThunderboltOutlined,
  ExperimentOutlined,
  ReloadOutlined,
} from '@ant-design/icons'
import type { QCReport, QCMetric, QCStatus } from '@/types/ai'
import { aiService } from '@/services/ai.service'
import { DesignTokens } from '@/styles/design-tokens'
import './AutoQCPanel.css'

const { Text, Paragraph, Title } = Typography
const { Panel } = Collapse

interface AutoQCPanelProps {
  sampleId: string
  sampleName?: string
  fastqFiles?: string[]
  bamFile?: string
  onRerun?: () => void
}

export const AutoQCPanel: React.FC<AutoQCPanelProps> = ({
  sampleId,
  sampleName,
  fastqFiles,
  bamFile,
  onRerun,
}) => {
  const [loading, setLoading] = useState(false)
  const [report, setReport] = useState<QCReport | null>(null)
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    runQCAnalysis()
  }, [sampleId])

  const runQCAnalysis = async () => {
    setLoading(true)
    setError(null)

    try {
      const result = await aiService.runAutoQC({
        sampleId,
        fastqFiles,
        bamFile,
        includeVisualizations: true,
      })
      setReport(result)
    } catch (err) {
      setError('Failed to run QC analysis')
      console.error('QC analysis error:', err)
    } finally {
      setLoading(false)
    }
  }

  const getStatusColor = (status: QCStatus): string => {
    const colors: Record<QCStatus, string> = {
      excellent: DesignTokens.colors.success.main,
      good: '#52c41a',
      acceptable: DesignTokens.colors.warning.main,
      poor: '#ff7a45',
      failed: DesignTokens.colors.error.main,
    }
    return colors[status] || DesignTokens.colors.gray[500]
  }

  const getStatusIcon = (status: QCStatus) => {
    if (status === 'excellent' || status === 'good') {
      return <CheckCircleOutlined style={{ color: getStatusColor(status) }} />
    }
    if (status === 'acceptable') {
      return <WarningOutlined style={{ color: getStatusColor(status) }} />
    }
    return <CloseCircleOutlined style={{ color: getStatusColor(status) }} />
  }

  const getSeverityColor = (severity: 'critical' | 'warning' | 'info'): string => {
    const colors = {
      critical: DesignTokens.colors.error.main,
      warning: DesignTokens.colors.warning.main,
      info: DesignTokens.colors.info.main,
    }
    return colors[severity]
  }

  if (loading) {
    return (
      <Card className="auto-qc-card">
        <div style={{ textAlign: 'center', padding: '48px' }}>
          <Spin size="large" />
          <Paragraph type="secondary" style={{ marginTop: 16 }}>
            AI is analyzing sample quality...
          </Paragraph>
          <Text type="secondary" style={{ fontSize: 12 }}>
            This may take a few moments
          </Text>
        </div>
      </Card>
    )
  }

  if (error) {
    return (
      <Card className="auto-qc-card">
        <Alert
          message="QC Analysis Error"
          description={error}
          type="error"
          showIcon
          action={
            <Button size="small" onClick={runQCAnalysis}>
              <ReloadOutlined /> Retry
            </Button>
          }
        />
      </Card>
    )
  }

  if (!report) {
    return null
  }

  return (
    <Card
      className="auto-qc-card"
      title={
        <Space>
          <ExperimentOutlined style={{ color: DesignTokens.colors.primary.main }} />
          <Text strong>AI Quality Control Report</Text>
          <Tag color="blue">
            <ThunderboltOutlined /> Auto-generated
          </Tag>
        </Space>
      }
      extra={
        <Button
          size="small"
          icon={<ReloadOutlined />}
          onClick={() => {
            runQCAnalysis()
            onRerun?.()
          }}
        >
          Rerun
        </Button>
      }
    >
      {/* Overall Status */}
      <div className="qc-overview">
        <Row gutter={[16, 16]}>
          <Col xs={24} sm={12} md={6}>
            <Card className="qc-stat-card">
              <Statistic
                title="Overall Score"
                value={report.overallScore}
                suffix="/ 100"
                valueStyle={{ color: getStatusColor(report.overallStatus) }}
                prefix={getStatusIcon(report.overallStatus)}
              />
              <Tag
                color={getStatusColor(report.overallStatus)}
                style={{ marginTop: 8, width: '100%', textAlign: 'center' }}
              >
                {report.overallStatus.toUpperCase()}
              </Tag>
            </Card>
          </Col>

          <Col xs={24} sm={12} md={6}>
            <Card className="qc-stat-card">
              <Statistic
                title="Passed Checks"
                value={report.passedChecks}
                suffix={`/ ${report.totalChecks}`}
                valueStyle={{
                  color:
                    report.passedChecks === report.totalChecks
                      ? DesignTokens.colors.success.main
                      : DesignTokens.colors.warning.main,
                }}
                prefix={<CheckCircleOutlined />}
              />
              <Progress
                percent={(report.passedChecks / report.totalChecks) * 100}
                strokeColor={
                  report.passedChecks === report.totalChecks
                    ? DesignTokens.colors.success.main
                    : DesignTokens.colors.warning.main
                }
                showInfo={false}
                style={{ marginTop: 8 }}
              />
            </Card>
          </Col>

          <Col xs={24} sm={12} md={6}>
            <Card className="qc-stat-card">
              <Statistic
                title="Sample"
                value={sampleName || report.sampleName}
                valueStyle={{ fontSize: 16 }}
              />
              <Text type="secondary" style={{ fontSize: 12 }}>
                ID: {report.sampleId}
              </Text>
            </Card>
          </Col>

          <Col xs={24} sm={12} md={6}>
            <Card className="qc-stat-card">
              <Statistic
                title="Analysis Time"
                value={new Date(report.timestamp).toLocaleString()}
                valueStyle={{ fontSize: 14 }}
              />
            </Card>
          </Col>
        </Row>
      </div>

      {/* Issues */}
      {report.issues && report.issues.length > 0 && (
        <div className="qc-issues">
          <Title level={5}>
            <WarningOutlined /> Issues Detected
          </Title>
          <Space direction="vertical" style={{ width: '100%' }} size="small">
            {report.issues.map((issue, idx) => (
              <Alert
                key={idx}
                message={issue.message}
                description={issue.suggestion}
                type={issue.severity === 'critical' ? 'error' : issue.severity}
                showIcon
                icon={
                  <span style={{ color: getSeverityColor(issue.severity) }}>
                    {issue.severity === 'critical' ? (
                      <CloseCircleOutlined />
                    ) : issue.severity === 'warning' ? (
                      <WarningOutlined />
                    ) : (
                      <InfoCircleOutlined />
                    )}
                  </span>
                }
              />
            ))}
          </Space>
        </div>
      )}

      {/* QC Metrics */}
      <div className="qc-metrics">
        <Title level={5}>Quality Metrics</Title>
        <Collapse defaultActiveKey={report.metrics.filter((m) => m.status !== 'excellent').map((m) => m.name)}>
          {report.metrics.map((metric) => (
            <Panel
              key={metric.name}
              header={
                <Space style={{ width: '100%', justifyContent: 'space-between' }}>
                  <Space>
                    {getStatusIcon(metric.status)}
                    <Text strong>{metric.name}</Text>
                  </Space>
                  <Space>
                    <Text code>
                      {metric.value} {metric.unit || ''}
                    </Text>
                    <Tag color={getStatusColor(metric.status)}>
                      {metric.status.toUpperCase()}
                    </Tag>
                  </Space>
                </Space>
              }
            >
              <div className="metric-details">
                <Paragraph>{metric.explanation}</Paragraph>

                {metric.threshold && (
                  <div className="metric-thresholds">
                    <Text strong style={{ display: 'block', marginBottom: 8 }}>
                      Thresholds:
                    </Text>
                    <Space direction="vertical" size="small" style={{ width: '100%' }}>
                      <div className="threshold-item">
                        <Tag color="green">Excellent</Tag>
                        <Text>≥ {metric.threshold.excellent}</Text>
                      </div>
                      <div className="threshold-item">
                        <Tag color="lime">Good</Tag>
                        <Text>≥ {metric.threshold.good}</Text>
                      </div>
                      <div className="threshold-item">
                        <Tag color="orange">Acceptable</Tag>
                        <Text>≥ {metric.threshold.acceptable}</Text>
                      </div>
                    </Space>
                  </div>
                )}

                <Progress
                  percent={
                    metric.threshold
                      ? Math.min((metric.value / metric.threshold.excellent) * 100, 100)
                      : 0
                  }
                  strokeColor={getStatusColor(metric.status)}
                  style={{ marginTop: 12 }}
                />
              </div>
            </Panel>
          ))}
        </Collapse>
      </div>

      {/* Recommendations */}
      {report.recommendations && report.recommendations.length > 0 && (
        <div className="qc-recommendations">
          <Alert
            message={
              <Space>
                <InfoCircleOutlined />
                <Text strong>AI Recommendations</Text>
              </Space>
            }
            description={
              <ul style={{ margin: 0, paddingLeft: 20 }}>
                {report.recommendations.map((rec, idx) => (
                  <li key={idx}>{rec}</li>
                ))}
              </ul>
            }
            type="info"
            showIcon={false}
          />
        </div>
      )}
    </Card>
  )
}
