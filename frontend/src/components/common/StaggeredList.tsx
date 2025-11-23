/**
 * StaggeredList - 交错动画列表
 * 列表项按顺序淡入，无需外部动画库
 */
import type React from 'react'
import { Children } from 'react'
import { FadeIn } from './FadeIn'

export interface StaggeredListProps {
  children: React.ReactNode
  /**
   * 每项之间的延迟 (ms)
   */
  staggerDelay?: number
  /**
   * 基础延迟 (ms)
   */
  baseDelay?: number
  /**
   * 动画持续时间 (ms)
   */
  duration?: number
  /**
   * 动画方向
   */
  direction?: 'up' | 'down' | 'left' | 'right' | 'none'
  className?: string
}

export const StaggeredList: React.FC<StaggeredListProps> = ({
  children,
  staggerDelay = 50,
  baseDelay = 0,
  duration = 250,
  direction = 'up',
  className = '',
}) => {
  const childArray = Children.toArray(children)

  return (
    <div className={className}>
      {childArray.map((child, index) => (
        <FadeIn
          key={index}
          delay={baseDelay + index * staggerDelay}
          duration={duration}
          direction={direction}
          once={true}
        >
          {child}
        </FadeIn>
      ))}
    </div>
  )
}

export default StaggeredList
