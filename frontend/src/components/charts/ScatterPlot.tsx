/**
 * Scatter Plot Component
 */
import React from 'react'
import { Chart } from './Chart'
import type { EChartsOption } from 'echarts'
import { designTokens } from '@/config'

export interface ScatterPlotProps {
  data: {
    x: number[]
    y: number[]
    name?: string
    labels?: string[]
  }[]
  title?: string
  xAxisLabel?: string
  yAxisLabel?: string
  height?: number
  showRegression?: boolean
  loading?: boolean
}

export const ScatterPlot: React.FC<ScatterPlotProps> = ({
  data,
  title,
  xAxisLabel,
  yAxisLabel,
  height = 400,
  showRegression = false,
  loading = false,
}) => {
  const option: EChartsOption = {
    tooltip: {
      trigger: 'item',
      formatter: (params: any) => {
        const { seriesName, data, dataIndex } = params
        const label = data[2] || `Point ${dataIndex + 1}`
        return `${seriesName}<br/>${label}<br/>X: ${data[0]}<br/>Y: ${data[1]}`
      },
    },
    legend: {
      data: data.map((d) => d.name || 'Series'),
      top: 10,
    },
    grid: {
      left: '3%',
      right: '7%',
      bottom: '10%',
      containLabel: true,
    },
    xAxis: {
      type: 'value',
      name: xAxisLabel,
      nameLocation: 'middle',
      nameGap: 30,
      splitLine: {
        lineStyle: {
          type: 'dashed',
        },
      },
    },
    yAxis: {
      type: 'value',
      name: yAxisLabel,
      nameLocation: 'middle',
      nameGap: 50,
      splitLine: {
        lineStyle: {
          type: 'dashed',
        },
      },
    },
    series: data.map((item, index) => ({
      name: item.name || `Series ${index + 1}`,
      type: 'scatter',
      data: item.x.map((xVal, i) => [xVal, item.y[i], item.labels?.[i]]),
      symbolSize: 8,
      itemStyle: {
        color: [
          designTokens.colors.primary.default,
          designTokens.colors.success.default,
          designTokens.colors.warning.default,
          designTokens.colors.error.default,
        ][index % 4],
        opacity: 0.7,
      },
      emphasis: {
        itemStyle: {
          opacity: 1,
          shadowBlur: 10,
          shadowOffsetX: 0,
          shadowOffsetY: 0,
          shadowColor: 'rgba(0, 0, 0, 0.5)',
        },
      },
    })),
  }

  return <Chart option={option} title={title} height={height} loading={loading} />
}

export default ScatterPlot
