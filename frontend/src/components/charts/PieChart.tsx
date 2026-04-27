/**
 * Pie / Donut Chart — Plotly implementation.
 */
import type React from 'react'
import { useMemo } from 'react'
import { Chart } from './Chart'
import type { Data, Layout } from 'plotly.js'

export interface PieChartProps {
  data: {
    name: string
    value: number
  }[]
  title?: string
  height?: number
  /**
   * Inner radius (0-1). Pass a value > 0 for a donut chart. Defaults to
   * 0.4 to roughly match the previous ECharts radius=['40%', '70%'].
   */
  hole?: number
  /** Render as nightingale rose (variable radius). */
  roseType?: boolean
  loading?: boolean
}

export const PieChart: React.FC<PieChartProps> = ({
  data,
  title,
  height = 400,
  hole = 0.4,
  roseType = false,
  loading = false,
}) => {
  const traces = useMemo<Data[]>(
    () => [
      {
        type: 'pie',
        labels: data.map((d) => d.name),
        values: data.map((d) => d.value),
        hole: roseType ? 0 : hole,
        textinfo: 'label+percent',
        textposition: 'auto',
        hovertemplate: '%{label}<br>%{value} (%{percent})<extra></extra>',
      } as Data,
    ],
    [data, hole, roseType],
  )

  const layout = useMemo<Partial<Layout>>(() => ({ showlegend: true }), [])

  return <Chart data={traces} layout={layout} title={title} height={height} loading={loading} />
}

export default PieChart
