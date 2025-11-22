/**
 * Design Tokens for NGSmodule
 * Unified design system values accessible in TypeScript/JavaScript
 *
 * These tokens mirror the CSS variables in variables.css
 * Use these when you need to access design values in JS/TS code
 *
 * @example
 * import { DesignTokens } from '@/styles/design-tokens'
 *
 * const cardStyle = {
 *   borderRadius: DesignTokens.radius.lg,
 *   padding: DesignTokens.spacing.lg,
 *   boxShadow: DesignTokens.shadow.md
 * }
 */

export const DesignTokens = {
  /**
   * Brand Colors
   * Primary color palette for the platform
   */
  colors: {
    // Primary - Scientific Blue
    primary: {
      main: '#2563eb',
      hover: '#1d4ed8',
      active: '#1e40af',
      light: '#dbeafe',
      lighter: '#eff6ff',
    },

    // Semantic Colors
    success: {
      main: '#10b981',
      hover: '#059669',
      light: '#d1fae5',
    },

    warning: {
      main: '#f59e0b',
      hover: '#d97706',
      light: '#fef3c7',
    },

    error: {
      main: '#ef4444',
      hover: '#dc2626',
      light: '#fee2e2',
    },

    info: {
      main: '#06b6d4',
      hover: '#0891b2',
      light: '#cffafe',
    },

    // Neutral Colors
    white: '#ffffff',
    black: '#000000',

    gray: {
      50: '#f9fafb',
      100: '#f3f4f6',
      200: '#e5e7eb',
      300: '#d1d5db',
      400: '#9ca3af',
      500: '#6b7280',
      600: '#4b5563',
      700: '#374151',
      800: '#1f2937',
      900: '#111827',
    },

    // Contextual Colors
    background: {
      primary: '#ffffff',
      secondary: '#f9fafb',
      tertiary: '#f3f4f6',
      hover: '#f3f4f6',
      active: '#e5e7eb',
      overlay: 'rgba(0, 0, 0, 0.45)',
    },

    text: {
      primary: '#111827',
      secondary: '#4b5563',
      tertiary: '#6b7280',
      disabled: '#9ca3af',
      inverse: '#ffffff',
    },

    border: {
      default: '#e5e7eb',
      dark: '#d1d5db',
      light: '#f3f4f6',
    },
  },

  /**
   * Spacing Scale
   * Consistent spacing values for layout and components
   */
  spacing: {
    xs: 4,
    sm: 8,
    md: 16,
    lg: 24,
    xl: 32,
    '2xl': 48,
    '3xl': 64,
  },

  /**
   * Typography
   * Font families, sizes, weights, and line heights
   */
  typography: {
    fontFamily: {
      sans: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
      mono: '"Fira Code", "Courier New", monospace',
    },

    fontSize: {
      xs: 12,
      sm: 14,
      md: 16,
      lg: 18,
      xl: 20,
      '2xl': 24,
      '3xl': 30,
      '4xl': 36,
    },

    fontWeight: {
      normal: 400,
      medium: 500,
      semibold: 600,
      bold: 700,
    },

    lineHeight: {
      tight: 1.25,
      normal: 1.5,
      relaxed: 1.75,
    },
  },

  /**
   * Border Radius
   * Rounding values for UI elements
   */
  radius: {
    sm: 2,
    md: 4,
    lg: 8,
    xl: 12,
    '2xl': 16,
    full: 9999,
  },

  /**
   * Shadows
   * Elevation system for depth and hierarchy
   */
  shadow: {
    sm: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
    md: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
    lg: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
    xl: '0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
    '2xl': '0 25px 50px -12px rgba(0, 0, 0, 0.25)',
    none: 'none',
  },

  /**
   * Z-Index Layers
   * Stacking order for overlays and modals
   */
  zIndex: {
    dropdown: 1000,
    sticky: 1020,
    fixed: 1030,
    modalBackdrop: 1040,
    modal: 1050,
    popover: 1060,
    tooltip: 1070,
  },

  /**
   * Transitions
   * Animation timing functions
   */
  transition: {
    fast: '150ms cubic-bezier(0.4, 0, 0.2, 1)',
    base: '250ms cubic-bezier(0.4, 0, 0.2, 1)',
    slow: '350ms cubic-bezier(0.4, 0, 0.2, 1)',
  },

  /**
   * Animation Durations
   * Timing values for animations
   */
  duration: {
    fast: 150,
    base: 250,
    slow: 350,
  },

  /**
   * Card & Container
   * Default styling for card components
   */
  card: {
    background: '#ffffff',
    shadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
    radius: 8,
    padding: 24,
  },

  /**
   * Layout Dimensions
   * Fixed sizes for layout components
   */
  layout: {
    sidebarWidth: 240,
    sidebarCollapsedWidth: 64,
    headerHeight: 64,
    footerHeight: 48,
  },

  /**
   * Breakpoints
   * Responsive design breakpoints (in pixels)
   */
  breakpoints: {
    xs: 480,
    sm: 576,
    md: 768,
    lg: 992,
    xl: 1200,
    '2xl': 1600,
  },
} as const

/**
 * Type-safe design tokens
 * Extracts types from the design tokens object
 */
export type DesignTokensType = typeof DesignTokens

/**
 * Helper function to generate px values
 * @param value - Numeric value to convert to px
 */
export const px = (value: number): string => `${value}px`

/**
 * Helper function to get spacing value in px
 * @param size - Spacing size key
 */
export const spacing = (size: keyof typeof DesignTokens.spacing): string => {
  return px(DesignTokens.spacing[size])
}

/**
 * Helper function to get font size in px
 * @param size - Font size key
 */
export const fontSize = (size: keyof typeof DesignTokens.typography.fontSize): string => {
  return px(DesignTokens.typography.fontSize[size])
}

/**
 * Helper function to get border radius in px
 * @param size - Radius size key
 */
export const borderRadius = (size: keyof typeof DesignTokens.radius): string | number => {
  const value = DesignTokens.radius[size]
  return value === 9999 ? value : px(value)
}

/**
 * Dark mode color overrides
 * Use these when implementing dark mode
 */
export const DarkModeTokens = {
  colors: {
    background: {
      primary: '#111827',
      secondary: '#1f2937',
      tertiary: '#374151',
    },

    text: {
      primary: '#f3f4f6',
      secondary: '#d1d5db',
      tertiary: '#9ca3af',
    },

    border: {
      default: '#374151',
    },

    card: {
      background: '#1f2937',
    },
  },
} as const

/**
 * Animation Easings
 * Common easing functions for animations
 */
export const Easings = {
  easeIn: 'cubic-bezier(0.4, 0, 1, 1)',
  easeOut: 'cubic-bezier(0, 0, 0.2, 1)',
  easeInOut: 'cubic-bezier(0.4, 0, 0.2, 1)',
  linear: 'linear',
  bounce: 'cubic-bezier(0.68, -0.55, 0.265, 1.55)',
} as const

/**
 * Common gradient patterns
 */
export const Gradients = {
  primary: `linear-gradient(135deg, ${DesignTokens.colors.primary.main} 0%, ${DesignTokens.colors.primary.active} 100%)`,
  primaryVertical: `linear-gradient(180deg, ${DesignTokens.colors.primary.main} 0%, ${DesignTokens.colors.primary.active} 100%)`,
  success: `linear-gradient(135deg, ${DesignTokens.colors.success.main} 0%, ${DesignTokens.colors.success.hover} 100%)`,
  warning: `linear-gradient(135deg, ${DesignTokens.colors.warning.main} 0%, ${DesignTokens.colors.warning.hover} 100%)`,
  error: `linear-gradient(135deg, ${DesignTokens.colors.error.main} 0%, ${DesignTokens.colors.error.hover} 100%)`,
  glass: 'linear-gradient(135deg, rgba(255, 255, 255, 0.1) 0%, rgba(255, 255, 255, 0.05) 100%)',
} as const

export default DesignTokens
