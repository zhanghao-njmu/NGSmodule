import { useState } from 'react'
import { Card, Row, Col, Typography, Space, Tabs, Tag, Statistic, Alert } from 'antd'
import {
  ThunderboltOutlined,
  BulbOutlined,
  ExperimentOutlined,
  WarningOutlined,
  GroupOutlined,
  RobotOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
} from '@ant-design/icons'
import { PageHeader, StatisticCard, EnhancedEmptyState, FadeIn, StaggeredList } from '@/components/common'
import { DesignTokens } from '@/styles/design-tokens'
import './AIDashboard.css'

const { Title, Paragraph } = Typography
const { TabPane } = Tabs

export const AIDashboard: React.FC = () => {
  const [activeTab, setActiveTab] = useState('overview')

  // Mock statistics - replace with real data
  const stats = {
    recommendationsGenerated: 245,
    qcReportsGenerated: 189,
    anomaliesDetected: 12,
    groupingsSuggested: 34,
    timesSaved: 127, // hours
    accuracyRate: 94.2,
  }

  return (
    <div className="ai-dashboard">
      {/* Header with animation */}
      <FadeIn direction="up" delay={0} duration={300}>
        <PageHeader
          title="AI Intelligence Center"
          subtitle="AI-powered insights and automation for your NGS workflow"
          icon={<RobotOutlined />}
        />
      </FadeIn>

      {/* Stats Overview with staggered animation */}
      <StaggeredList staggerDelay={60} baseDelay={100} direction="up">
        <Row gutter={[16, 16]} style={{ marginBottom: 24 }}>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="Recommendations"
              value={stats.recommendationsGenerated}
              icon={<BulbOutlined />}
              color={DesignTokens.colors.warning.main}
            />
          </Col>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="QC Reports"
              value={stats.qcReportsGenerated}
              icon={<ExperimentOutlined />}
              color={DesignTokens.colors.primary.main}
            />
          </Col>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="Anomalies Found"
              value={stats.anomaliesDetected}
              icon={<WarningOutlined />}
              color={DesignTokens.colors.error.main}
            />
          </Col>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="Smart Groupings"
              value={stats.groupingsSuggested}
              icon={<GroupOutlined />}
              color={DesignTokens.colors.success.main}
            />
          </Col>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="Time Saved"
              value={stats.timesSaved}
              suffix="hrs"
              icon={<ClockCircleOutlined />}
              color={DesignTokens.colors.info.main}
            />
          </Col>
          <Col xs={12} sm={8} md={4}>
            <StatisticCard
              title="Accuracy Rate"
              value={stats.accuracyRate}
              suffix="%"
              icon={<CheckCircleOutlined />}
              color={DesignTokens.colors.success.main}
            />
          </Col>
        </Row>
      </StaggeredList>

      {/* AI Features Tabs with animation */}
      <FadeIn direction="up" delay={150} duration={300}>
        <Card className="ai-features-card">
          <Tabs activeKey={activeTab} onChange={setActiveTab} size="large">
            {/* Overview */}
            <TabPane
              tab={
                <Space>
                  <ThunderboltOutlined />
                  Overview
                </Space>
              }
              key="overview"
            >
              <div className="ai-overview">
                <Row gutter={[24, 24]}>
                  {/* Parameter Recommendations */}
                  <Col xs={24} md={12}>
                    <Card className="feature-card hover-lift" hoverable onClick={() => setActiveTab('recommendations')}>
                      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
                        <div className="feature-icon">
                          <BulbOutlined style={{ fontSize: 48, color: DesignTokens.colors.warning.main }} />
                        </div>
                        <div>
                          <Title level={4}>Parameter Recommendations</Title>
                          <Paragraph type="secondary">
                            AI-powered suggestions for optimal pipeline parameters based on your data characteristics,
                            similar successful runs, and best practices.
                          </Paragraph>
                        </div>
                        <div className="feature-stats">
                          <Statistic
                            title="Recommendations Generated"
                            value={stats.recommendationsGenerated}
                            prefix={<BulbOutlined />}
                          />
                        </div>
                        <Tag color="blue">
                          <ThunderboltOutlined /> Auto-optimized
                        </Tag>
                      </Space>
                    </Card>
                  </Col>

                  {/* Auto QC */}
                  <Col xs={24} md={12}>
                    <Card className="feature-card hover-lift" hoverable onClick={() => setActiveTab('qc')}>
                      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
                        <div className="feature-icon">
                          <ExperimentOutlined style={{ fontSize: 48, color: DesignTokens.colors.primary.main }} />
                        </div>
                        <div>
                          <Title level={4}>Auto Quality Control</Title>
                          <Paragraph type="secondary">
                            Automated quality assessment with intelligent scoring, issue detection, and actionable
                            recommendations for sample quality improvement.
                          </Paragraph>
                        </div>
                        <div className="feature-stats">
                          <Statistic
                            title="QC Reports Generated"
                            value={stats.qcReportsGenerated}
                            prefix={<ExperimentOutlined />}
                          />
                        </div>
                        <Tag color="green">
                          <CheckCircleOutlined /> Automated
                        </Tag>
                      </Space>
                    </Card>
                  </Col>

                  {/* Anomaly Detection */}
                  <Col xs={24} md={12}>
                    <Card className="feature-card hover-lift" hoverable onClick={() => setActiveTab('anomalies')}>
                      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
                        <div className="feature-icon">
                          <WarningOutlined style={{ fontSize: 48, color: DesignTokens.colors.error.main }} />
                        </div>
                        <div>
                          <Title level={4}>Anomaly Detection</Title>
                          <Paragraph type="secondary">
                            Intelligent detection of quality drops, contamination, batch effects, and outliers with
                            automatic alerts and fix suggestions.
                          </Paragraph>
                        </div>
                        <div className="feature-stats">
                          <Statistic
                            title="Anomalies Detected"
                            value={stats.anomaliesDetected}
                            prefix={<WarningOutlined />}
                          />
                        </div>
                        <Tag color="red">
                          <WarningOutlined /> Real-time
                        </Tag>
                      </Space>
                    </Card>
                  </Col>

                  {/* Smart Grouping */}
                  <Col xs={24} md={12}>
                    <Card className="feature-card hover-lift" hoverable onClick={() => setActiveTab('grouping')}>
                      <Space direction="vertical" size="middle" style={{ width: '100%' }}>
                        <div className="feature-icon">
                          <GroupOutlined style={{ fontSize: 48, color: DesignTokens.colors.success.main }} />
                        </div>
                        <div>
                          <Title level={4}>Smart Sample Grouping</Title>
                          <Paragraph type="secondary">
                            Automatic sample grouping based on metadata, QC metrics, and expression patterns with
                            suggested comparisons.
                          </Paragraph>
                        </div>
                        <div className="feature-stats">
                          <Statistic
                            title="Groupings Suggested"
                            value={stats.groupingsSuggested}
                            prefix={<GroupOutlined />}
                          />
                        </div>
                        <Tag color="cyan">
                          <ThunderboltOutlined /> ML-powered
                        </Tag>
                      </Space>
                    </Card>
                  </Col>
                </Row>

                {/* Benefits Section */}
                <Card style={{ marginTop: 24 }} className="benefits-card">
                  <Title level={4}>
                    <ThunderboltOutlined /> AI Intelligence Benefits
                  </Title>
                  <Row gutter={[16, 16]}>
                    <Col xs={24} md={8}>
                      <Space direction="vertical">
                        <CheckCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.success.main }} />
                        <Title level={5}>Time Savings</Title>
                        <Paragraph type="secondary">
                          Save hours of manual parameter tuning and quality checking with automated AI recommendations.
                        </Paragraph>
                      </Space>
                    </Col>
                    <Col xs={24} md={8}>
                      <Space direction="vertical">
                        <CheckCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.success.main }} />
                        <Title level={5}>Higher Accuracy</Title>
                        <Paragraph type="secondary">
                          Leverage best practices and historical data for {stats.accuracyRate}% accuracy in
                          recommendations.
                        </Paragraph>
                      </Space>
                    </Col>
                    <Col xs={24} md={8}>
                      <Space direction="vertical">
                        <CheckCircleOutlined style={{ fontSize: 24, color: DesignTokens.colors.success.main }} />
                        <Title level={5}>Proactive Detection</Title>
                        <Paragraph type="secondary">
                          Catch issues before they become problems with real-time anomaly detection and alerts.
                        </Paragraph>
                      </Space>
                    </Col>
                  </Row>
                </Card>
              </div>
            </TabPane>

            {/* Parameter Recommendations Tab */}
            <TabPane
              tab={
                <Space>
                  <BulbOutlined />
                  Recommendations
                </Space>
              }
              key="recommendations"
            >
              <EnhancedEmptyState
                type="noData"
                title="No recommendations yet"
                description="Parameter recommendations appear when creating or editing pipelines"
                action={{
                  text: 'Go to Pipelines',
                  onClick: () => {
                    /* navigate to pipelines */
                  },
                  icon: <ThunderboltOutlined />,
                }}
                size="default"
              />
            </TabPane>

            {/* QC Tab */}
            <TabPane
              tab={
                <Space>
                  <ExperimentOutlined />
                  Quality Control
                </Space>
              }
              key="qc"
            >
              <EnhancedEmptyState
                type="noData"
                title="No QC reports yet"
                description="QC reports are generated automatically for each sample"
                action={{
                  text: 'View Samples',
                  onClick: () => {
                    /* navigate to samples */
                  },
                  icon: <ExperimentOutlined />,
                }}
                size="default"
              />
            </TabPane>

            {/* Anomalies Tab */}
            <TabPane
              tab={
                <Space>
                  <WarningOutlined />
                  Anomalies
                  {stats.anomaliesDetected > 0 && <Tag color="red">{stats.anomaliesDetected}</Tag>}
                </Space>
              }
              key="anomalies"
            >
              <Alert
                message="Anomaly Detection Active"
                description="Real-time monitoring is enabled for all active items. You will be notified of any detected anomalies."
                type="info"
                showIcon
                style={{ marginBottom: 16 }}
              />
              <EnhancedEmptyState
                type="noData"
                title="No anomalies detected"
                description="No anomalies detected in your recent data"
                size="default"
              />
            </TabPane>

            {/* Grouping Tab */}
            <TabPane
              tab={
                <Space>
                  <GroupOutlined />
                  Sample Grouping
                </Space>
              }
              key="grouping"
            >
              <EnhancedEmptyState
                type="noData"
                title="No groupings yet"
                description="Smart grouping suggestions are available in project analysis"
                action={{
                  text: 'Go to Projects',
                  onClick: () => {
                    /* navigate to items */
                  },
                  icon: <GroupOutlined />,
                }}
                size="default"
              />
            </TabPane>
          </Tabs>
        </Card>
      </FadeIn>
    </div>
  )
}
