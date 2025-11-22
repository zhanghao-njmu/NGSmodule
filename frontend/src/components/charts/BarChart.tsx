/**
 * Bar Chart Component
 */
import React from 'react'
import { Chart } from './Chart'
import type { EChartsOption } from 'echarts'
import { designTokens } from '@/config'

export interface BarChartProps {
  data: {
    categories: string[]
    values: number[]
    name?: string
  }[]
  title?: string
  xAxisLabel?: string
  yAxisLabel?: string
  height?: number
  horizontal?: boolean
  stacked?: boolean
  loading?: boolean
}

export const BarChart: React.FC<BarChartProps> = ({
  data,
  title,
  xAxisLabel,
  yAxisLabel,
  height = 400,
  horizontal = false,
  stacked = false,
  loading = false,
}) => {
  const option: EChartsOption = {
    tooltip: {
      trigger: 'axis',
      axisPointer: {
        type: 'shadow',
      },
    },
    legend: {
      data: data.map((d) => d.name || 'Series'),
      top: 10,
    },
    grid: {
      left: '3%',
      right: '4%',
      bottom: '10%',
      containLabel: true,
    },
    xAxis: {
      type: horizontal ? 'value' : 'category',
      data: horizontal ? undefined : data[0]?.categories || [],
      name: xAxisLabel,
      nameLocation: 'middle',
      nameGap: 30,
    },
    yAxis: {
      type: horizontal ? 'category' : 'value',
      data: horizontal ? data[0]?.categories || [] : undefined,
      name: yAxisLabel,
      nameLocation: 'middle',
      nameGap: 50,
    },
    series: data.map((item, index) => ({
      name: item.name || `Series ${index + 1}`,
      type: 'bar',
      data: item.values,
      stack: stacked ? 'total' : undefined,
      itemStyle: {
        color: [
          designTokens.colors.primary.default,
          designTokens.colors.success.default,
          designTokens.colors.warning.default,
          designTokens.colors.error.default,
        ][index % 4],
      },
    })),
  }

  return <Chart option={option} title={title} height={height} loading={loading} />
}

export default BarChart
