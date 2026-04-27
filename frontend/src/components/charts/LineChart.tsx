/**
 * Line Chart — Plotly implementation.
 * External API kept identical to the previous ECharts version.
 */
import React, { useMemo } from 'react'
import { Chart } from './Chart'
import type { Data, Layout } from 'plotly.js'
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

const PALETTE = [
  designTokens.colors.primary.default,
  designTokens.colors.success.default,
  designTokens.colors.warning.default,
  designTokens.colors.error.default,
]

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
  const traces = useMemo<Data[]>(
    () =>
      data.map((series, index) => ({
        type: 'scatter',
        mode: 'lines+markers',
        x: series.x as any,
        y: series.y,
        name: series.name || `Series ${index + 1}`,
        line: {
          shape: smooth ? 'spline' : 'linear',
          color: PALETTE[index % PALETTE.length],
          width: 2,
        },
        ...(showArea ? { fill: 'tozeroy', fillcolor: `${PALETTE[index % PALETTE.length]}33` } : {}),
      })),
    [data, smooth, showArea],
  )

  const layout = useMemo<Partial<Layout>>(
    () => ({
      xaxis: { title: xAxisLabel, automargin: true },
      yaxis: { title: yAxisLabel, automargin: true },
    }),
    [xAxisLabel, yAxisLabel],
  )

  return <Chart data={traces} layout={layout} title={title} height={height} loading={loading} />
}

export default LineChart
