/**
 * Line Chart Component
 */
import React from 'react'
import { Chart } from './Chart'
import type { EChartsOption } from 'echarts'
import { designTokens } from '@/config'

export interface LineChartProps {
  data: {
    x: (string | number)[]
    y: number[]
    name?: string
  }[]
  title?: string
  xAxisLabel?: string
  yAxisLabel?: string
  height?: number
  smooth?: boolean
  showArea?: boolean
  loading?: boolean
}

export const LineChart: React.FC<LineChartProps> = ({
  data,
  title,
  xAxisLabel,
  yAxisLabel,
  height = 400,
  smooth = true,
  showArea = false,
  loading = false,
}) => {
  const option: EChartsOption = {
    tooltip: {
      trigger: 'axis',
      axisPointer: {
        type: 'cross',
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
      type: 'category',
      boundaryGap: false,
      data: data[0]?.x || [],
      name: xAxisLabel,
      nameLocation: 'middle',
      nameGap: 30,
    },
    yAxis: {
      type: 'value',
      name: yAxisLabel,
      nameLocation: 'middle',
      nameGap: 50,
    },
    series: data.map((item, index) => ({
      name: item.name || `Series ${index + 1}`,
      type: 'line',
      smooth,
      data: item.y,
      areaStyle: showArea ? {} : undefined,
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

export default LineChart
