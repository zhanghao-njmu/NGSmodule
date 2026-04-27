/**
 * Generic Chart wrapper around react-plotly.js.
 *
 * Plotly is the de-facto standard for scientific visualization (volcano
 * plots, Manhattan plots, heatmaps), so we standardize on it across the
 * app and drop ECharts. The wrapper exposes a small typed surface so
 * downstream chart components don't import plotly types directly.
 */
import React, { useMemo } from 'react'
import Plot from 'react-plotly.js'
import type { Data, Layout, Config } from 'plotly.js'
import { Card, Spin } from 'antd'

export interface ChartProps {
  data: Data[]
  layout?: Partial<Layout>
  config?: Partial<Config>
  title?: string
  height?: number | string
  loading?: boolean
  className?: string
  style?: React.CSSProperties
}

const DEFAULT_CONFIG: Partial<Config> = {
  responsive: true,
  displaylogo: false,
  modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
}

const DEFAULT_LAYOUT: Partial<Layout> = {
  autosize: true,
  margin: { l: 50, r: 30, t: 40, b: 50 },
  font: {
    family:
      "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
    size: 12,
  },
  paper_bgcolor: 'transparent',
  plot_bgcolor: 'transparent',
  legend: {
    orientation: 'h',
    yanchor: 'bottom',
    y: 1.02,
    xanchor: 'right',
    x: 1,
  },
  hovermode: 'closest',
}

export const Chart: React.FC<ChartProps> = ({
  data,
  layout = {},
  config = {},
  title,
  height = 400,
  loading = false,
  className,
  style,
}) => {
  const mergedLayout = useMemo<Partial<Layout>>(
    () => ({
      ...DEFAULT_LAYOUT,
      ...layout,
      autosize: true,
    }),
    [layout],
  )

  const mergedConfig = useMemo<Partial<Config>>(() => ({ ...DEFAULT_CONFIG, ...config }), [config])

  const heightCss = typeof height === 'number' ? `${height}px` : height

  const content = (
    <Spin spinning={loading}>
      <div className={className} style={{ width: '100%', height: heightCss, ...style }}>
        <Plot
          data={data}
          layout={mergedLayout}
          config={mergedConfig}
          useResizeHandler
          style={{ width: '100%', height: '100%' }}
        />
      </div>
    </Spin>
  )

  if (title) {
    return (
      <Card title={title} bordered={false}>
        {content}
      </Card>
    )
  }

  return content
}

export default Chart
