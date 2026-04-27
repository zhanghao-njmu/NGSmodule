/**
 * Animation Utilities
 * 动画工具函数和配置
 */
import { designTokens } from '@/config/design-tokens'

/**
 * 动画持续时间
 */
export const animationDurations = {
  fast: designTokens.durations.fast,
  base: designTokens.durations.base,
  slow: designTokens.durations.slow,
} as const

/**
 * 动画缓动函数
 */
export const easings = {
  easeInOut: 'cubic-bezier(0.4, 0, 0.2, 1)',
  easeOut: 'cubic-bezier(0.0, 0, 0.2, 1)',
  easeIn: 'cubic-bezier(0.4, 0, 1, 1)',
  sharp: 'cubic-bezier(0.4, 0, 0.6, 1)',
} as const

/**
 * 页面过渡动画变体
 */
export const pageTransitions = {
  // 淡入淡出
  fade: {
    initial: { opacity: 0 },
    animate: { opacity: 1 },
    exit: { opacity: 0 },
    transition: { duration: animationDurations.base / 1000 },
  },

  // 从下向上滑入
  slideUp: {
    initial: { opacity: 0, y: 20 },
    animate: { opacity: 1, y: 0 },
    exit: { opacity: 0, y: -20 },
    transition: { duration: animationDurations.base / 1000, ease: easings.easeInOut },
  },

  // 从右向左滑入
  slideLeft: {
    initial: { opacity: 0, x: 20 },
    animate: { opacity: 1, x: 0 },
    exit: { opacity: 0, x: -20 },
    transition: { duration: animationDurations.base / 1000, ease: easings.easeInOut },
  },

  // 缩放
  scale: {
    initial: { opacity: 0, scale: 0.95 },
    animate: { opacity: 1, scale: 1 },
    exit: { opacity: 0, scale: 0.95 },
    transition: { duration: animationDurations.base / 1000, ease: easings.easeOut },
  },
} as const

/**
 * 列表项交错动画
 */
export const staggerChildren = {
  animate: {
    transition: {
      staggerChildren: 0.05,
    },
  },
}

export const listItemVariants = {
  initial: { opacity: 0, y: 10 },
  animate: { opacity: 1, y: 0 },
  exit: { opacity: 0, y: -10 },
}

/**
 * Modal 动画变体
 */
export const modalVariants = {
  initial: { opacity: 0, scale: 0.9, y: -20 },
  animate: { opacity: 1, scale: 1, y: 0 },
  exit: { opacity: 0, scale: 0.9, y: -20 },
  transition: { duration: animationDurations.base / 1000, ease: easings.easeOut },
}

/**
 * Hover 缩放动画
 */
export const hoverScale = {
  whileHover: { scale: 1.02 },
  whileTap: { scale: 0.98 },
  transition: { duration: animationDurations.fast / 1000 },
}

/**
 * CSS 过渡类名
 */
export const transitionClasses = {
  fade: {
    enter: 'transition-opacity duration-300 ease-in-out',
    enterFrom: 'opacity-0',
    enterTo: 'opacity-100',
    leave: 'transition-opacity duration-200 ease-in-out',
    leaveFrom: 'opacity-100',
    leaveTo: 'opacity-0',
  },
  slideUp: {
    enter: 'transition-all duration-300 ease-out',
    enterFrom: 'opacity-0 translate-y-4',
    enterTo: 'opacity-100 translate-y-0',
    leave: 'transition-all duration-200 ease-in',
    leaveFrom: 'opacity-100 translate-y-0',
    leaveTo: 'opacity-0 -translate-y-4',
  },
  scale: {
    enter: 'transition-all duration-300 ease-out',
    enterFrom: 'opacity-0 scale-95',
    enterTo: 'opacity-100 scale-100',
    leave: 'transition-all duration-200 ease-in',
    leaveFrom: 'opacity-100 scale-100',
    leaveTo: 'opacity-0 scale-95',
  },
}

/**
 * 创建CSS过渡样式
 */
export const createTransition = (
  properties: string[],
  duration: keyof typeof animationDurations = 'base',
  easing: keyof typeof easings = 'easeInOut',
): string => {
  const durationMs = animationDurations[duration]
  const easingFn = easings[easing]
  return properties.map((prop) => `${prop} ${durationMs}ms ${easingFn}`).join(', ')
}

/**
 * 延迟执行
 */
export const delay = (ms: number): Promise<void> => {
  return new Promise((resolve) => setTimeout(resolve, ms))
}

/**
 * 生成交错延迟
 */
export const getStaggerDelay = (index: number, baseDelay: number = 50): number => {
  return index * baseDelay
}
