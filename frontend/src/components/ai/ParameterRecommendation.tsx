import { useState, useEffect } from 'react'
import { Card, Space, Typography, Tag, Button, Tooltip, Progress, Alert, Spin } from 'antd'
import {
  BulbOutlined,
  CheckCircleOutlined,
  InfoCircleOutlined,
  ThunderboltOutlined,
  LikeOutlined,
  DislikeOutlined,
} from '@ant-design/icons'
import type { ParameterRecommendation, PipelineRecommendation } from '@/types/ai'
import { aiService } from '@/services/ai.service'
import { DesignTokens } from '@/styles/design-tokens'
import './ParameterRecommendation.css'

const { Text, Paragraph, Title } = Typography

interface ParameterRecommendationProps {
  pipelineType: string
  sampleType?: string
  organism?: string
  onApplyRecommendation?: (parameters: Record<string, any>) => void
  onParameterChange?: (parameter: string, value: any) => void
}

export const ParameterRecommendationWidget: React.FC<ParameterRecommendationProps> = ({
  pipelineType,
  sampleType,
  organism,
  onApplyRecommendation,
  onParameterChange,
}) => {
  const [loading, setLoading] = useState(false)
  const [recommendation, setRecommendation] = useState<PipelineRecommendation | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [expandedParams, setExpandedParams] = useState<Set<string>>(new Set())

  useEffect(() => {
    fetchRecommendations()
  }, [pipelineType, sampleType, organism])

  const fetchRecommendations = async () => {
    if (!pipelineType) return

    setLoading(true)
    setError(null)

    try {
      const result = await aiService.getParameterRecommendations({
        pipelineType,
        sampleType,
        organism,
      })
      setRecommendation(result)
    } catch (err) {
      setError('Failed to fetch AI recommendations')
      console.error('Recommendation error:', err)
    } finally {
      setLoading(false)
    }
  }

  const handleApplyAll = () => {
    if (!recommendation) return

    const parameters: Record<string, any> = {}
    recommendation.parameters.forEach((param) => {
      parameters[param.parameter] = param.recommendedValue
    })

    onApplyRecommendation?.(parameters)
  }

  const handleApplyParameter = (param: ParameterRecommendation) => {
    onParameterChange?.(param.parameter, param.recommendedValue)
  }

  const toggleParameterDetails = (paramName: string) => {
    const newExpanded = new Set(expandedParams)
    if (newExpanded.has(paramName)) {
      newExpanded.delete(paramName)
    } else {
      newExpanded.add(paramName)
    }
    setExpandedParams(newExpanded)
  }

  const getConfidenceColor = (confidence: number): string => {
    if (confidence >= 0.8) return DesignTokens.colors.success.main
    if (confidence >= 0.6) return DesignTokens.colors.warning.main
    return DesignTokens.colors.error.main
  }

  const getConfidenceLabel = (confidence: number): string => {
    if (confidence >= 0.8) return 'High Confidence'
    if (confidence >= 0.6) return 'Medium Confidence'
    return 'Low Confidence'
  }

  if (loading) {
    return (
      <Card className="ai-recommendation-card">
        <div style={{ textAlign: 'center', padding: '24px' }}>
          <Spin size="large" />
          <Paragraph type="secondary" style={{ marginTop: 16 }}>
            AI is analyzing optimal parameters...
          </Paragraph>
        </div>
      </Card>
    )
  }

  if (error) {
    return (
      <Card className="ai-recommendation-card">
        <Alert
          message="Recommendation Error"
          description={error}
          type="error"
          showIcon
          action={
            <Button size="small" onClick={fetchRecommendations}>
              Retry
            </Button>
          }
        />
      </Card>
    )
  }

  if (!recommendation) {
    return null
  }

  return (
    <Card
      className="ai-recommendation-card"
      title={
        <Space>
          <BulbOutlined style={{ color: DesignTokens.colors.warning.main }} />
          <Text strong>AI-Powered Parameter Recommendations</Text>
        </Space>
      }
      extra={
        <Space>
          <Tag color="blue">
            <ThunderboltOutlined /> AI Optimized
          </Tag>
          <Button type="primary" size="small" onClick={handleApplyAll}>
            Apply All
          </Button>
        </Space>
      }
    >
      {/* Overall Confidence */}
      <div className="recommendation-overview">
        <Space direction="vertical" style={{ width: '100%' }} size="small">
          <div className="confidence-bar">
            <Text type="secondary">Overall Confidence</Text>
            <Progress
              percent={recommendation.overallConfidence * 100}
              strokeColor={getConfidenceColor(recommendation.overallConfidence)}
              format={(percent) => `${percent?.toFixed(0)}%`}
            />
          </div>

          {recommendation.estimatedRuntime && (
            <Text type="secondary">
              <InfoCircleOutlined /> Estimated Runtime: {recommendation.estimatedRuntime}
            </Text>
          )}
        </Space>
      </div>

      {/* Warnings */}
      {recommendation.warnings && recommendation.warnings.length > 0 && (
        <Alert
          message="Warnings"
          description={
            <ul style={{ margin: 0, paddingLeft: 20 }}>
              {recommendation.warnings.map((warning, idx) => (
                <li key={idx}>{warning}</li>
              ))}
            </ul>
          }
          type="warning"
          showIcon
          style={{ marginTop: 16 }}
        />
      )}

      {/* Parameter Recommendations */}
      <div className="parameter-recommendations">
        <Title level={5} style={{ marginTop: 16 }}>
          Recommended Parameters
        </Title>

        <Space direction="vertical" style={{ width: '100%' }} size="middle">
          {recommendation.parameters.map((param) => (
            <Card
              key={param.parameter}
              size="small"
              className="parameter-card"
              hoverable
              onClick={() => toggleParameterDetails(param.parameter)}
            >
              <div className="parameter-header">
                <Space style={{ flex: 1 }}>
                  <CheckCircleOutlined
                    style={{ color: getConfidenceColor(param.confidence) }}
                  />
                  <div>
                    <Text strong>{param.parameter}</Text>
                    <br />
                    <Text code>{String(param.recommendedValue)}</Text>
                  </div>
                </Space>

                <Space>
                  <Tooltip title={getConfidenceLabel(param.confidence)}>
                    <Tag color={getConfidenceColor(param.confidence)}>
                      {(param.confidence * 100).toFixed(0)}%
                    </Tag>
                  </Tooltip>
                  <Button
                    type="primary"
                    size="small"
                    onClick={(e) => {
                      e.stopPropagation()
                      handleApplyParameter(param)
                    }}
                  >
                    Apply
                  </Button>
                </Space>
              </div>

              {expandedParams.has(param.parameter) && (
                <div className="parameter-details">
                  <Paragraph type="secondary">{param.reason}</Paragraph>

                  {param.basedOn && (
                    <div className="based-on">
                      <Text type="secondary" style={{ fontSize: 12 }}>
                        Based on:{' '}
                      </Text>
                      {param.basedOn.similarProjects && (
                        <Tag>{param.basedOn.similarProjects} similar projects</Tag>
                      )}
                      {param.basedOn.bestPractices && <Tag>Best practices</Tag>}
                      {param.basedOn.publicDatasets && (
                        <Tag>{param.basedOn.publicDatasets} public datasets</Tag>
                      )}
                    </div>
                  )}

                  {param.alternatives && param.alternatives.length > 0 && (
                    <div className="alternatives">
                      <Text strong style={{ fontSize: 12 }}>
                        Alternatives:
                      </Text>
                      <Space direction="vertical" size="small" style={{ width: '100%' }}>
                        {param.alternatives.map((alt, idx) => (
                          <div key={idx} className="alternative-option">
                            <Space>
                              <Text code>{String(alt.value)}</Text>
                              <Tag size="small">
                                {(alt.confidence * 100).toFixed(0)}%
                              </Tag>
                            </Space>
                            <Text type="secondary" style={{ fontSize: 11 }}>
                              {alt.reason}
                            </Text>
                          </div>
                        ))}
                      </Space>
                    </div>
                  )}
                </div>
              )}
            </Card>
          ))}
        </Space>
      </div>

      {/* Suggestions */}
      {recommendation.suggestions && recommendation.suggestions.length > 0 && (
        <Alert
          message="Suggestions"
          description={
            <ul style={{ margin: 0, paddingLeft: 20 }}>
              {recommendation.suggestions.map((suggestion, idx) => (
                <li key={idx}>{suggestion}</li>
              ))}
            </ul>
          }
          type="info"
          showIcon
          style={{ marginTop: 16 }}
        />
      )}

      {/* Feedback */}
      <div className="recommendation-feedback">
        <Text type="secondary" style={{ fontSize: 12 }}>
          Was this recommendation helpful?
        </Text>
        <Space>
          <Button
            size="small"
            icon={<LikeOutlined />}
            onClick={() => {
              // TODO: Submit positive feedback
            }}
          >
            Yes
          </Button>
          <Button
            size="small"
            icon={<DislikeOutlined />}
            onClick={() => {
              // TODO: Submit negative feedback
            }}
          >
            No
          </Button>
        </Space>
      </div>
    </Card>
  )
}
