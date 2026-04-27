/**
 * Bar Chart — Plotly implementation.
 */
import type React from 'react'
import { useMemo } from 'react'
import { Chart } from './Chart'
import type { Data, Layout } from 'plotly.js'
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

const PALETTE = [
  designTokens.colors.primary.default,
  designTokens.colors.success.default,
  designTokens.colors.warning.default,
  designTokens.colors.error.default,
]

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
  const traces = useMemo<Data[]>(
    () =>
      data.map((series, index) => ({
        type: 'bar',
        orientation: horizontal ? 'h' : 'v',
        x: horizontal ? series.values : series.categories,
        y: horizontal ? series.categories : series.values,
        name: series.name || `Series ${index + 1}`,
        marker: { color: PALETTE[index % PALETTE.length] },
      })),
    [data, horizontal],
  )

  const layout = useMemo<Partial<Layout>>(
    () => ({
      barmode: stacked ? 'stack' : 'group',
      xaxis: { title: xAxisLabel, automargin: true },
      yaxis: { title: yAxisLabel, automargin: true },
    }),
    [xAxisLabel, yAxisLabel, stacked],
  )

  return <Chart data={traces} layout={layout} title={title} height={height} loading={loading} />
}

export default BarChart
