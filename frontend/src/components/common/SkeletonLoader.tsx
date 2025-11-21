/**
 * Skeleton Loader Component - Beautiful loading placeholders
 */
import React from 'react'
import { Card, Skeleton, Row, Col } from 'antd'

interface SkeletonLoaderProps {
  type?: 'table' | 'card' | 'statistics' | 'form'
  rows?: number
  cardCount?: number
}

export const SkeletonLoader: React.FC<SkeletonLoaderProps> = ({
  type = 'table',
  rows = 5,
  cardCount = 4,
}) => {
  if (type === 'statistics') {
    return (
      <Row gutter={16} style={{ marginBottom: 24 }}>
        {Array.from({ length: cardCount }).map((_, index) => (
          <Col xs={24} sm={12} lg={6} key={index}>
            <Card>
              <Skeleton.Input active style={{ width: '100%', height: 80 }} />
            </Card>
          </Col>
        ))}
      </Row>
    )
  }

  if (type === 'card') {
    return (
      <Row gutter={[16, 16]}>
        {Array.from({ length: cardCount }).map((_, index) => (
          <Col xs={24} sm={12} lg={8} key={index}>
            <Card>
              <Skeleton active paragraph={{ rows: 3 }} />
            </Card>
          </Col>
        ))}
      </Row>
    )
  }

  if (type === 'form') {
    return (
      <Card>
        <Skeleton active paragraph={{ rows: 6 }} />
      </Card>
    )
  }

  // Default: table skeleton
  return (
    <Card>
      <Skeleton active paragraph={{ rows }} />
    </Card>
  )
}

export default SkeletonLoader
