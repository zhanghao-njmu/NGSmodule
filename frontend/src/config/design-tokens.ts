/**
 * Design Tokens - TypeScript Configuration
 * 与 CSS 变量保持同步，用于 JS/TS 运行时访问
 */

export const designTokens = {
  // ========== Colors ==========
  colors: {
    // Brand Colors
    primary: {
      default: '#2563eb',
      hover: '#1d4ed8',
      active: '#1e40af',
      light: '#dbeafe',
      lighter: '#eff6ff',
    },
    success: {
      default: '#10b981',
      hover: '#059669',
      light: '#d1fae5',
    },
    warning: {
      default: '#f59e0b',
      hover: '#d97706',
      light: '#fef3c7',
    },
    error: {
      default: '#ef4444',
      hover: '#dc2626',
      light: '#fee2e2',
    },
    info: {
      default: '#06b6d4',
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

    // Semantic Colors
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

  // ========== Spacing ==========
  spacing: {
    xs: 4,
    sm: 8,
    md: 16,
    lg: 24,
    xl: 32,
    '2xl': 48,
    '3xl': 64,
  },

  // ========== Typography ==========
  typography: {
    fontFamily: {
      sans: "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif",
      mono: "'Fira Code', 'Courier New', monospace",
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

  // ========== Border Radius ==========
  borderRadius: {
    sm: 2,
    md: 4,
    lg: 8,
    xl: 12,
    '2xl': 16,
    full: 9999,
  },

  // ========== Shadows ==========
  shadows: {
    sm: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',
    md: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
    lg: '0 10px 15px -3px rgba(0, 0, 0, 0.1), 0 4px 6px -2px rgba(0, 0, 0, 0.05)',
    xl: '0 20px 25px -5px rgba(0, 0, 0, 0.1), 0 10px 10px -5px rgba(0, 0, 0, 0.04)',
    '2xl': '0 25px 50px -12px rgba(0, 0, 0, 0.25)',
  },

  // ========== Z-Index ==========
  zIndex: {
    dropdown: 1000,
    sticky: 1020,
    fixed: 1030,
    modalBackdrop: 1040,
    modal: 1050,
    popover: 1060,
    tooltip: 1070,
  },

  // ========== Transitions ==========
  transitions: {
    fast: '150ms cubic-bezier(0.4, 0, 0.2, 1)',
    base: '250ms cubic-bezier(0.4, 0, 0.2, 1)',
    slow: '350ms cubic-bezier(0.4, 0, 0.2, 1)',
  },

  // ========== Animation Durations ==========
  durations: {
    fast: 150,
    base: 250,
    slow: 350,
  },

  // ========== Layout ==========
  layout: {
    sidebarWidth: 240,
    sidebarCollapsedWidth: 64,
    headerHeight: 64,
    footerHeight: 48,
  },

  // ========== Breakpoints ==========
  breakpoints: {
    xs: 480,
    sm: 576,
    md: 768,
    lg: 992,
    xl: 1200,
    '2xl': 1600,
  },
} as const

// Type exports for auto-completion
export type DesignTokens = typeof designTokens
export type ColorKey = keyof typeof designTokens.colors
export type SpacingKey = keyof typeof designTokens.spacing

// Utility function to access tokens
export const getToken = <T extends keyof DesignTokens>(
  category: T
): DesignTokens[T] => {
  return designTokens[category]
}

// Helper to get CSS variable
export const getCSSVar = (varName: string): string => {
  return `var(--${varName})`
}

// Status color mapping
export const statusColors = {
  pending: designTokens.colors.gray[400],
  running: designTokens.colors.primary.default,
  completed: designTokens.colors.success.default,
  failed: designTokens.colors.error.default,
  paused: designTokens.colors.warning.default,
  cancelled: designTokens.colors.gray[500],
} as const

// File type color mapping
export const fileTypeColors = {
  fastq: designTokens.colors.primary.default,
  bam: designTokens.colors.success.default,
  vcf: '#722ed1', // purple
  sam: designTokens.colors.warning.default,
  bed: designTokens.colors.info.default,
  gff: '#eb2f96', // magenta
} as const
