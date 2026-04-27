/**
 * Scatter Plot — Plotly implementation.
 */
import type React from 'react'
import { useMemo } from 'react'
import { Chart } from './Chart'
import type { Data, Layout } from 'plotly.js'
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

const PALETTE = [
  designTokens.colors.primary.default,
  designTokens.colors.success.default,
  designTokens.colors.warning.default,
  designTokens.colors.error.default,
]

/** Fit y = m*x + b via simple least squares. */
function fitLinearRegression(xs: number[], ys: number[]): { slope: number; intercept: number } {
  const n = xs.length
  if (n === 0) {
    return { slope: 0, intercept: 0 }
  }
  let sumX = 0
  let sumY = 0
  let sumXY = 0
  let sumXX = 0
  for (let i = 0; i < n; i++) {
    sumX += xs[i]
    sumY += ys[i]
    sumXY += xs[i] * ys[i]
    sumXX += xs[i] * xs[i]
  }
  const denom = n * sumXX - sumX * sumX
  const slope = denom === 0 ? 0 : (n * sumXY - sumX * sumY) / denom
  const intercept = (sumY - slope * sumX) / n
  return { slope, intercept }
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
  const traces = useMemo<Data[]>(() => {
    const out: Data[] = []
    data.forEach((series, index) => {
      const color = PALETTE[index % PALETTE.length]
      out.push({
        type: 'scatter',
        mode: 'markers',
        x: series.x,
        y: series.y,
        name: series.name || `Series ${index + 1}`,
        text: series.labels,
        marker: { color, size: 8, opacity: 0.7 },
        hovertemplate: '%{text}<br>x=%{x}, y=%{y}<extra></extra>',
      })

      if (showRegression && series.x.length >= 2) {
        const { slope, intercept } = fitLinearRegression(series.x, series.y)
        const xMin = Math.min(...series.x)
        const xMax = Math.max(...series.x)
        out.push({
          type: 'scatter',
          mode: 'lines',
          x: [xMin, xMax],
          y: [slope * xMin + intercept, slope * xMax + intercept],
          name: `${series.name || `Series ${index + 1}`} fit`,
          line: { color, dash: 'dash', width: 2 },
          showlegend: false,
          hoverinfo: 'skip',
        })
      }
    })
    return out
  }, [data, showRegression])

  const layout = useMemo<Partial<Layout>>(
    () => ({
      xaxis: { title: xAxisLabel, automargin: true, zeroline: false },
      yaxis: { title: yAxisLabel, automargin: true, zeroline: false },
    }),
    [xAxisLabel, yAxisLabel],
  )

  return <Chart data={traces} layout={layout} title={title} height={height} loading={loading} />
}

export default ScatterPlot
