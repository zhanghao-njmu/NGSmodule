/**
 * Generic Chart Component - Wrapper for ECharts
 */
import React, { useRef, useEffect } from 'react'
import * as echarts from 'echarts'
import { Card } from 'antd'
import type { EChartsOption } from 'echarts'

export interface ChartProps {
  option: EChartsOption
  title?: string
  height?: number | string
  loading?: boolean
  className?: string
  style?: React.CSSProperties
  onChartReady?: (chart: echarts.ECharts) => void
}

export const Chart: React.FC<ChartProps> = ({
  option,
  title,
  height = 400,
  loading = false,
  className,
  style,
  onChartReady,
}) => {
  const chartRef = useRef<HTMLDivElement>(null)
  const chartInstance = useRef<echarts.ECharts | null>(null)

  useEffect(() => {
    if (!chartRef.current) return

    // Initialize chart
    if (!chartInstance.current) {
      chartInstance.current = echarts.init(chartRef.current)
      onChartReady?.(chartInstance.current)
    }

    // Set option
    chartInstance.current.setOption(option, true)

    // Handle loading
    if (loading) {
      chartInstance.current.showLoading()
    } else {
      chartInstance.current.hideLoading()
    }

    // Handle resize
    const handleResize = () => {
      chartInstance.current?.resize()
    }
    window.addEventListener('resize', handleResize)

    return () => {
      window.removeEventListener('resize', handleResize)
    }
  }, [option, loading, onChartReady])

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      chartInstance.current?.dispose()
    }
  }, [])

  const content = (
    <div
      ref={chartRef}
      className={className}
      style={{
        width: '100%',
        height: typeof height === 'number' ? `${height}px` : height,
        ...style,
      }}
    />
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
