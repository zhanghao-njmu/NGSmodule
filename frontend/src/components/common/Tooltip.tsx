/**
 * Enhanced Tooltip Component - Beautiful tooltips with animations
 */
import React from 'react'
import { Tooltip as AntTooltip } from 'antd'
import type { TooltipProps as AntTooltipProps } from 'antd'

interface TooltipProps extends AntTooltipProps {
  variant?: 'default' | 'info' | 'success' | 'warning' | 'error'
}

export const Tooltip: React.FC<TooltipProps> = ({
  variant = 'default',
  overlayInnerStyle,
  ...props
}) => {
  const variantColors = {
    default: '#1f2937',
    info: '#2563eb',
    success: '#10b981',
    warning: '#f59e0b',
    error: '#ef4444',
  }

  return (
    <AntTooltip
      mouseEnterDelay={0.3}
      overlayInnerStyle={{
        backgroundColor: variantColors[variant],
        borderRadius: 6,
        padding: '6px 12px',
        fontSize: 13,
        boxShadow: '0 4px 12px rgba(0, 0, 0, 0.15)',
        ...overlayInnerStyle,
      }}
      {...props}
    />
  )
}

export default Tooltip
