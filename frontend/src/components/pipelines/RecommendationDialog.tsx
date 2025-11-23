/**
 * AI Parameter Recommendation Dialog
 * Visualizes parameter recommendations with confidence scores
 */
import type React from 'react'
import {
  Modal,
  Card,
  Row,
  Col,
  Statistic,
  Progress,
  Tag,
  Descriptions,
  Alert,
  Space,
  Typography,
  Divider,
  Timeline,
  Button,
} from 'antd'
import {
  BulbOutlined,
  ThunderboltOutlined,
  CheckCircleOutlined,
  HistoryOutlined,
  RocketOutlined,
  InfoCircleOutlined,
} from '@ant-design/icons'
import type { ParameterRecommendationResponse } from '@/types/pipeline'

const { Title, Text, Paragraph } = Typography

interface RecommendationDialogProps {
  open: boolean
  recommendation: ParameterRecommendationResponse | null
  loading: boolean
  onClose: () => void
  onApply: (params: Record<string, any>) => void
}

export const RecommendationDialog: React.FC<RecommendationDialogProps> = ({
  open,
  recommendation,
  loading,
  onClose,
  onApply,
}) => {
  if (!recommendation && !loading) {
    return null
  }

  // Calculate confidence level
  const getConfidenceLevel = (score: number) => {
    if (score >= 0.8) {
      return { level: 'High', color: 'success', status: 'active' }
    }
    if (score >= 0.6) {
      return { level: 'Medium', color: 'normal', status: 'active' }
    }
    return { level: 'Low', color: 'exception', status: 'exception' }
  }

  const confidenceConfig = recommendation
    ? getConfidenceLevel(recommendation.confidence_score)
    : { level: 'Unknown', color: 'normal', status: 'active' }

  // Group parameters by type
  const groupParametersByType = (params: Record<string, any>) => {
    const groups: Record<string, any[]> = {
      numeric: [],
      boolean: [],
      string: [],
    }

    Object.entries(params).forEach(([key, value]) => {
      if (typeof value === 'number') {
        groups.numeric.push({ key, value, type: 'number' })
      } else if (typeof value === 'boolean') {
        groups.boolean.push({ key, value, type: 'boolean' })
      } else {
        groups.string.push({ key, value, type: 'string' })
      }
    })

    return groups
  }

  const paramGroups = recommendation
    ? groupParametersByType(recommendation.recommended_params)
    : { numeric: [], boolean: [], string: [] }

  const handleApply = () => {
    if (recommendation) {
      onApply(recommendation.recommended_params)
      onClose()
    }
  }

  return (
    <Modal
      title={
        <Space>
          <BulbOutlined style={{ color: '#faad14', fontSize: 20 }} />
          <Title level={4} style={{ margin: 0 }}>
            AI Parameter Recommendations
          </Title>
        </Space>
      }
      open={open}
      onCancel={onClose}
      width={800}
      footer={
        !loading && recommendation ? (
          <Space>
            <Button onClick={onClose}>Cancel</Button>
            <Button type="primary" icon={<ThunderboltOutlined />} onClick={handleApply}>
              Apply Recommendations
            </Button>
          </Space>
        ) : null
      }
    >
      {loading ? (
        <Card loading style={{ minHeight: 300 }}>
          <div style={{ textAlign: 'center', padding: '60px 0' }}>
            <BulbOutlined style={{ fontSize: 64, color: '#faad14', marginBottom: 16 }} />
            <Title level={4}>Analyzing Historical Tasks...</Title>
            <Paragraph type="secondary">
              Our AI is analyzing your previous successful pipeline runs to recommend optimal parameters.
            </Paragraph>
          </div>
        </Card>
      ) : recommendation ? (
        <Space direction="vertical" size="large" style={{ width: '100%' }}>
          {/* Confidence Score */}
          <Row gutter={16}>
            <Col span={8}>
              <Card size="small">
                <Statistic
                  title="Confidence Score"
                  value={Math.round(recommendation.confidence_score * 100)}
                  suffix="%"
                  prefix={<CheckCircleOutlined />}
                  valueStyle={{
                    color:
                      confidenceConfig.color === 'success'
                        ? '#52c41a'
                        : confidenceConfig.color === 'normal'
                          ? '#1890ff'
                          : '#faad14',
                  }}
                />
                <Progress
                  percent={Math.round(recommendation.confidence_score * 100)}
                  size="small"
                  status={confidenceConfig.status as any}
                  showInfo={false}
                  style={{ marginTop: 8 }}
                />
                <Tag
                  color={
                    confidenceConfig.level === 'High'
                      ? 'success'
                      : confidenceConfig.level === 'Medium'
                        ? 'processing'
                        : 'warning'
                  }
                  style={{ marginTop: 8 }}
                >
                  {confidenceConfig.level} Confidence
                </Tag>
              </Card>
            </Col>

            <Col span={8}>
              <Card size="small">
                <Statistic
                  title="Based on Tasks"
                  value={recommendation.based_on_tasks}
                  prefix={<HistoryOutlined />}
                  valueStyle={{ color: '#722ed1' }}
                />
                <Text type="secondary" style={{ fontSize: 12 }}>
                  Analyzed successful runs
                </Text>
              </Card>
            </Col>

            <Col span={8}>
              <Card size="small">
                <Statistic
                  title="Parameters"
                  value={Object.keys(recommendation.recommended_params).length}
                  prefix={<RocketOutlined />}
                  valueStyle={{ color: '#13c2c2' }}
                />
                <Text type="secondary" style={{ fontSize: 12 }}>
                  Recommended settings
                </Text>
              </Card>
            </Col>
          </Row>

          {/* Explanation */}
          <Alert
            message="Recommendation Explanation"
            description={recommendation.explanation}
            type="info"
            showIcon
            icon={<InfoCircleOutlined />}
          />

          {/* Parameters */}
          <Card title="Recommended Parameters" size="small">
            <Timeline
              items={[
                ...(paramGroups.numeric.length > 0
                  ? [
                      {
                        dot: <ThunderboltOutlined style={{ color: '#1890ff' }} />,
                        children: (
                          <div>
                            <Text strong>Numeric Parameters</Text>
                            <Descriptions column={2} size="small" bordered style={{ marginTop: 8 }}>
                              {paramGroups.numeric.map(({ key, value }) => (
                                <Descriptions.Item label={key} key={key}>
                                  <Tag color="blue">{value}</Tag>
                                </Descriptions.Item>
                              ))}
                            </Descriptions>
                          </div>
                        ),
                      },
                    ]
                  : []),
                ...(paramGroups.boolean.length > 0
                  ? [
                      {
                        dot: <CheckCircleOutlined style={{ color: '#52c41a' }} />,
                        children: (
                          <div>
                            <Text strong>Boolean Flags</Text>
                            <Descriptions column={2} size="small" bordered style={{ marginTop: 8 }}>
                              {paramGroups.boolean.map(({ key, value }) => (
                                <Descriptions.Item label={key} key={key}>
                                  <Tag color={value ? 'success' : 'default'}>{value ? 'Enabled' : 'Disabled'}</Tag>
                                </Descriptions.Item>
                              ))}
                            </Descriptions>
                          </div>
                        ),
                      },
                    ]
                  : []),
                ...(paramGroups.string.length > 0
                  ? [
                      {
                        dot: <InfoCircleOutlined style={{ color: '#722ed1' }} />,
                        children: (
                          <div>
                            <Text strong>Text Parameters</Text>
                            <Descriptions column={1} size="small" bordered style={{ marginTop: 8 }}>
                              {paramGroups.string.map(({ key, value }) => (
                                <Descriptions.Item label={key} key={key}>
                                  <Tag color="purple">{value}</Tag>
                                </Descriptions.Item>
                              ))}
                            </Descriptions>
                          </div>
                        ),
                      },
                    ]
                  : []),
              ]}
            />
          </Card>

          <Divider />

          <Alert
            message="How it works"
            description={
              <ul style={{ marginBottom: 0, paddingLeft: 20 }}>
                <li>AI analyzes your previous successful pipeline runs</li>
                <li>Identifies parameter patterns that led to success</li>
                <li>Recommends the most frequently used values</li>
                <li>Confidence score based on data consistency and sample size</li>
              </ul>
            }
            type="success"
            showIcon
          />
        </Space>
      ) : null}
    </Modal>
  )
}

export default RecommendationDialog
