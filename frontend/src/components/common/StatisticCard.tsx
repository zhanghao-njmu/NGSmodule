/**
 * Statistic Card Component - Consistent statistics display
 */
import React from 'react'
import { Card, Statistic, Row, Col } from 'antd'
import type { StatisticProps } from 'antd'

export interface StatisticItem extends StatisticProps {
  key: string
  colSpan?: {
    xs?: number
    sm?: number
    md?: number
    lg?: number
    xl?: number
    xxl?: number
  }
}

interface StatisticCardProps {
  items: StatisticItem[]
  gutter?: number | [number, number]
  style?: React.CSSProperties
}

export const StatisticCard: React.FC<StatisticCardProps> = ({
  items,
  gutter = 16,
  style,
}) => {
  return (
    <Row gutter={gutter} style={{ marginBottom: 24, ...style }}>
      {items.map((item) => {
        const { key, colSpan = { xs: 24, sm: 12, lg: 6 }, ...statisticProps } = item
        return (
          <Col key={key} {...colSpan}>
            <Card>
              <Statistic {...statisticProps} />
            </Card>
          </Col>
        )
      })}
    </Row>
  )
}

export default StatisticCard
