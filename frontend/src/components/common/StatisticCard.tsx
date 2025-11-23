/**
 * Statistic Card Component - Consistent statistics display
 */
import type React from 'react'
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

export interface StatisticCardProps extends StatisticProps {
  // Single card mode props
  title?: string | React.ReactNode
  value?: string | number
  icon?: React.ReactNode
  color?: string
  suffix?: string
  prefix?: string

  // Multiple cards mode props
  items?: StatisticItem[]
  gutter?: number | [number, number]
  style?: React.CSSProperties
}

export const StatisticCard: React.FC<StatisticCardProps> = ({
  items,
  gutter = 16,
  style,
  title,
  value,
  icon,
  color,
  suffix,
  prefix,
  ...rest
}) => {
  // If items array is provided, render multiple statistics
  if (items) {
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

  // Otherwise, render a single statistic card
  return (
    <Card style={style} styles={{ body: { padding: '20px 24px' } }}>
      <Statistic title={title} value={value} suffix={suffix} prefix={icon || prefix} valueStyle={{ color }} {...rest} />
    </Card>
  )
}

export default StatisticCard
