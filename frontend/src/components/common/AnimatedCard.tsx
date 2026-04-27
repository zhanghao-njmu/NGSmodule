/**
 * Animated Card Component - Card with hover animations and effects
 */
import type React from 'react'
import { useState } from 'react'
import { Card } from 'antd'
import type { CardProps } from 'antd'

interface AnimatedCardProps extends CardProps {
  hoverEffect?: 'lift' | 'glow' | 'scale' | 'none'
  animationDelay?: number
}

export const AnimatedCard: React.FC<AnimatedCardProps> = ({
  hoverEffect = 'lift',
  animationDelay = 0,
  children,
  style,
  ...props
}) => {
  const [isHovered, setIsHovered] = useState(false)

  const getHoverStyle = () => {
    if (!isHovered || hoverEffect === 'none') {
      return {}
    }

    switch (hoverEffect) {
      case 'lift':
        return {
          transform: 'translateY(-4px)',
          boxShadow: '0 12px 24px -8px rgba(0, 0, 0, 0.15)',
        }
      case 'glow':
        return {
          boxShadow: '0 0 20px rgba(37, 99, 235, 0.3)',
        }
      case 'scale':
        return {
          transform: 'scale(1.02)',
        }
      default:
        return {}
    }
  }

  return (
    <Card
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
      style={{
        transition: 'all 0.3s cubic-bezier(0.4, 0, 0.2, 1)',
        animation: `fadeIn 0.5s ease-out ${animationDelay}ms both`,
        ...getHoverStyle(),
        ...style,
      }}
      {...props}
    >
      {children}
    </Card>
  )
}

export default AnimatedCard
