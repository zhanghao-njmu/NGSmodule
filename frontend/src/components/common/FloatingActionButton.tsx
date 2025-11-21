/**
 * Floating Action Button - Animated FAB for primary actions
 */
import React, { useState } from 'react'
import { Button, Tooltip } from 'antd'
import type { ButtonProps } from 'antd'
import { PlusOutlined } from '@ant-design/icons'

interface FloatingActionButtonProps extends ButtonProps {
  tooltip?: string
  position?: 'bottom-right' | 'bottom-left'
}

export const FloatingActionButton: React.FC<FloatingActionButtonProps> = ({
  tooltip = 'Add New',
  position = 'bottom-right',
  icon = <PlusOutlined />,
  ...props
}) => {
  const [isHovered, setIsHovered] = useState(false)

  const positionStyles = {
    'bottom-right': { bottom: 24, right: 24 },
    'bottom-left': { bottom: 24, left: 24 },
  }

  return (
    <Tooltip title={tooltip} placement="left">
      <Button
        type="primary"
        shape="circle"
        size="large"
        icon={icon}
        onMouseEnter={() => setIsHovered(true)}
        onMouseLeave={() => setIsHovered(false)}
        style={{
          position: 'fixed',
          ...positionStyles[position],
          width: 56,
          height: 56,
          boxShadow: isHovered
            ? '0 8px 16px rgba(37, 99, 235, 0.4)'
            : '0 4px 8px rgba(37, 99, 235, 0.3)',
          transform: isHovered ? 'scale(1.1)' : 'scale(1)',
          transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
          zIndex: 1000,
        }}
        {...props}
      />
    </Tooltip>
  )
}

export default FloatingActionButton
