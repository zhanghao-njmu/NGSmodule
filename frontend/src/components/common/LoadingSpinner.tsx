/**
 * Loading Spinner Component - Consistent loading state display
 */
import React from 'react'
import { Spin } from 'antd'
import { LoadingOutlined } from '@ant-design/icons'

interface LoadingSpinnerProps {
  size?: 'small' | 'default' | 'large'
  tip?: string
  fullscreen?: boolean
  style?: React.CSSProperties
}

export const LoadingSpinner: React.FC<LoadingSpinnerProps> = ({
  size = 'large',
  tip = 'Loading...',
  fullscreen = false,
  style,
}) => {
  const icon = <LoadingOutlined style={{ fontSize: size === 'large' ? 48 : 24 }} spin />

  if (fullscreen) {
    return (
      <div
        style={{
          display: 'flex',
          justifyContent: 'center',
          alignItems: 'center',
          minHeight: '100vh',
          ...style,
        }}
      >
        <Spin indicator={icon} size={size} tip={tip} />
      </div>
    )
  }

  return (
    <div className="loading-container" style={style}>
      <Spin indicator={icon} size={size} tip={tip} />
    </div>
  )
}

export default LoadingSpinner
