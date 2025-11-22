/**
 * FadeIn - 淡入动画组件
 * 简单的淡入效果，无需外部动画库
 */
import React, { useEffect, useState, useRef } from 'react'
import { designTokens } from '@/config/design-tokens'

export interface FadeInProps {
  children: React.ReactNode
  delay?: number
  duration?: number
  /**
   * 是否只在首次挂载时触发动画
   */
  once?: boolean
  /**
   * 动画方向
   */
  direction?: 'up' | 'down' | 'left' | 'right' | 'none'
  /**
   * 移动距离 (px)
   */
  distance?: number
  className?: string
  style?: React.CSSProperties
}

export const FadeIn: React.FC<FadeInProps> = ({
  children,
  delay = 0,
  duration = designTokens.durations.base,
  once = true,
  direction = 'up',
  distance = 20,
  className = '',
  style = {},
}) => {
  const [isVisible, setIsVisible] = useState(false)
  const [hasAnimated, setHasAnimated] = useState(false)
  const ref = useRef<HTMLDivElement>(null)

  useEffect(() => {
    // 如果已经动画过且 once=true，不再执行
    if (once && hasAnimated) return

    const observer = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting) {
          setTimeout(() => {
            setIsVisible(true)
            setHasAnimated(true)
          }, delay)
        } else if (!once) {
          setIsVisible(false)
        }
      },
      {
        threshold: 0.1,
        rootMargin: '50px',
      }
    )

    if (ref.current) {
      observer.observe(ref.current)
    }

    return () => {
      if (ref.current) {
        observer.unobserve(ref.current)
      }
    }
  }, [delay, once, hasAnimated])

  const getTransform = () => {
    if (direction === 'none' || isVisible) return 'none'

    switch (direction) {
      case 'up':
        return `translateY(${distance}px)`
      case 'down':
        return `translateY(-${distance}px)`
      case 'left':
        return `translateX(${distance}px)`
      case 'right':
        return `translateX(-${distance}px)`
      default:
        return 'none'
    }
  }

  return (
    <div
      ref={ref}
      className={className}
      style={{
        opacity: isVisible ? 1 : 0,
        transform: getTransform(),
        transition: `opacity ${duration}ms cubic-bezier(0.4, 0, 0.2, 1), transform ${duration}ms cubic-bezier(0.4, 0, 0.2, 1)`,
        ...style,
      }}
    >
      {children}
    </div>
  )
}

export default FadeIn
