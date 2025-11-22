/**
 * Pie Chart Component
 */
import React from 'react'
import { Chart } from './Chart'
import type { EChartsOption } from 'echarts'

export interface PieChartProps {
  data: {
    name: string
    value: number
  }[]
  title?: string
  height?: number
  radius?: string | [string, string]
  roseType?: boolean
  loading?: boolean
}

export const PieChart: React.FC<PieChartProps> = ({
  data,
  title,
  height = 400,
  radius = ['40%', '70%'],
  roseType = false,
  loading = false,
}) => {
  const option: EChartsOption = {
    tooltip: {
      trigger: 'item',
      formatter: '{a} <br/>{b}: {c} ({d}%)',
    },
    legend: {
      orient: 'vertical',
      left: 'left',
      data: data.map((d) => d.name),
    },
    series: [
      {
        name: title || 'Distribution',
        type: 'pie',
        radius,
        roseType: roseType ? 'radius' : undefined,
        data,
        emphasis: {
          itemStyle: {
            shadowBlur: 10,
            shadowOffsetX: 0,
            shadowColor: 'rgba(0, 0, 0, 0.5)',
          },
        },
        label: {
          formatter: '{b}: {d}%',
        },
      },
    ],
  }

  return <Chart option={option} title={title} height={height} loading={loading} />
}

export default PieChart
