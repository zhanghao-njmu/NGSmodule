/**
 * useMediaQuery Hook
 * 响应式媒体查询 Hook
 */
import { useState, useEffect } from 'react'
import { designTokens } from '@/config/design-tokens'

export type Breakpoint = 'xs' | 'sm' | 'md' | 'lg' | 'xl' | '2xl'

/**
 * 断点映射
 */
const breakpoints: Record<Breakpoint, number> = {
  xs: designTokens.breakpoints.xs,
  sm: designTokens.breakpoints.sm,
  md: designTokens.breakpoints.md,
  lg: designTokens.breakpoints.lg,
  xl: designTokens.breakpoints.xl,
  '2xl': designTokens.breakpoints['2xl'],
}

/**
 * useMediaQuery
 * 监听媒体查询变化
 *
 * @param query - 媒体查询字符串或断点名称
 * @returns 是否匹配查询
 *
 * @example
 * // 使用媒体查询字符串
 * const isMobile = useMediaQuery('(max-width: 768px)')
 *
 * // 使用断点名称
 * const isDesktop = useMediaQuery('lg')
 */
export const useMediaQuery = (query: string | Breakpoint): boolean => {
  // 如果是断点名称，转换为媒体查询
  const mediaQuery = query in breakpoints ? `(min-width: ${breakpoints[query as Breakpoint]}px)` : query

  const [matches, setMatches] = useState(() => {
    if (typeof window === 'undefined') {
      return false
    }
    return window.matchMedia(mediaQuery).matches
  })

  useEffect(() => {
    if (typeof window === 'undefined') {
      return
    }

    const mediaQueryList = window.matchMedia(mediaQuery)

    const handleChange = (event: MediaQueryListEvent) => {
      setMatches(event.matches)
    }

    // 现代浏览器
    if (mediaQueryList.addEventListener) {
      mediaQueryList.addEventListener('change', handleChange)
      return () => mediaQueryList.removeEventListener('change', handleChange)
    } else {
      // 旧浏览器兼容
      mediaQueryList.addListener(handleChange)
      return () => mediaQueryList.removeListener(handleChange)
    }
  }, [mediaQuery])

  return matches
}

/**
 * useBreakpoint
 * 获取当前断点
 *
 * @returns 当前断点名称
 *
 * @example
 * const breakpoint = useBreakpoint()
 * // 返回: 'xs' | 'sm' | 'md' | 'lg' | 'xl' | '2xl'
 */
export const useBreakpoint = (): Breakpoint => {
  const [breakpoint, setBreakpoint] = useState<Breakpoint>('md')

  useEffect(() => {
    if (typeof window === 'undefined') {
      return
    }

    const updateBreakpoint = () => {
      const width = window.innerWidth

      if (width < breakpoints.xs) {
        setBreakpoint('xs')
      } else if (width < breakpoints.sm) {
        setBreakpoint('xs')
      } else if (width < breakpoints.md) {
        setBreakpoint('sm')
      } else if (width < breakpoints.lg) {
        setBreakpoint('md')
      } else if (width < breakpoints.xl) {
        setBreakpoint('lg')
      } else if (width < breakpoints['2xl']) {
        setBreakpoint('xl')
      } else {
        setBreakpoint('2xl')
      }
    }

    updateBreakpoint()
    window.addEventListener('resize', updateBreakpoint)

    return () => window.removeEventListener('resize', updateBreakpoint)
  }, [])

  return breakpoint
}

/**
 * useIsMobile
 * 判断是否为移动设备
 *
 * @returns 是否为移动设备 (< md)
 */
export const useIsMobile = (): boolean => {
  return useMediaQuery(`(max-width: ${breakpoints.md - 1}px)`)
}

/**
 * useIsTablet
 * 判断是否为平板设备
 *
 * @returns 是否为平板设备 (md - lg)
 */
export const useIsTablet = (): boolean => {
  const isAboveMd = useMediaQuery('md')
  const isBelowLg = useMediaQuery(`(max-width: ${breakpoints.lg - 1}px)`)
  return isAboveMd && isBelowLg
}

/**
 * useIsDesktop
 * 判断是否为桌面设备
 *
 * @returns 是否为桌面设备 (>= lg)
 */
export const useIsDesktop = (): boolean => {
  return useMediaQuery('lg')
}

export default useMediaQuery
