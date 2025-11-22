/**
 * Configuration Index
 * 统一导出所有配置
 */

// Theme configuration
export { default as theme } from './theme'
export { antdTheme, darkTheme } from './theme.config'

// Design tokens
export { designTokens, getToken, getCSSVar, statusColors, fileTypeColors } from './design-tokens'
export type { DesignTokens, ColorKey, SpacingKey } from './design-tokens'
