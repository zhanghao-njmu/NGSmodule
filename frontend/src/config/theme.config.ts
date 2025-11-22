/**
 * Ant Design Theme Configuration
 * 使用 Design Tokens 配置 Ant Design 5 主题
 */
import type { ThemeConfig } from 'antd'
import { designTokens } from './design-tokens'

export const antdTheme: ThemeConfig = {
  token: {
    // ========== Brand Colors ==========
    colorPrimary: designTokens.colors.primary.default,
    colorSuccess: designTokens.colors.success.default,
    colorWarning: designTokens.colors.warning.default,
    colorError: designTokens.colors.error.default,
    colorInfo: designTokens.colors.info.default,

    // ========== Neutral Colors ==========
    colorBgBase: designTokens.colors.white,
    colorTextBase: designTokens.colors.text.primary,

    // ========== Border ==========
    colorBorder: designTokens.colors.border.default,
    colorBorderSecondary: designTokens.colors.border.light,

    // ========== Typography ==========
    fontFamily: designTokens.typography.fontFamily.sans,
    fontSize: designTokens.typography.fontSize.sm,
    fontSizeHeading1: designTokens.typography.fontSize['4xl'],
    fontSizeHeading2: designTokens.typography.fontSize['3xl'],
    fontSizeHeading3: designTokens.typography.fontSize['2xl'],
    fontSizeHeading4: designTokens.typography.fontSize.xl,
    fontSizeHeading5: designTokens.typography.fontSize.lg,

    // ========== Border Radius ==========
    borderRadius: designTokens.borderRadius.lg,
    borderRadiusSM: designTokens.borderRadius.md,
    borderRadiusLG: designTokens.borderRadius.xl,

    // ========== Spacing ==========
    padding: designTokens.spacing.md,
    paddingSM: designTokens.spacing.sm,
    paddingLG: designTokens.spacing.lg,
    paddingXL: designTokens.spacing.xl,

    margin: designTokens.spacing.md,
    marginSM: designTokens.spacing.sm,
    marginLG: designTokens.spacing.lg,
    marginXL: designTokens.spacing.xl,

    // ========== Line Height ==========
    lineHeight: designTokens.typography.lineHeight.normal,
    lineHeightHeading: designTokens.typography.lineHeight.tight,

    // ========== Control Heights ==========
    controlHeight: 36,
    controlHeightSM: 28,
    controlHeightLG: 44,

    // ========== Motion ==========
    motionDurationFast: `${designTokens.durations.fast}ms`,
    motionDurationMid: `${designTokens.durations.base}ms`,
    motionDurationSlow: `${designTokens.durations.slow}ms`,

    // ========== Z-Index ==========
    zIndexPopupBase: designTokens.zIndex.dropdown,
    zIndexBase: designTokens.zIndex.fixed,
  },

  components: {
    // ========== Button ==========
    Button: {
      borderRadius: designTokens.borderRadius.lg,
      controlHeight: 36,
      controlHeightLG: 44,
      controlHeightSM: 28,
      paddingContentHorizontal: designTokens.spacing.md,
      fontWeight: designTokens.typography.fontWeight.medium,
      primaryShadow: '0 2px 4px rgba(37, 99, 235, 0.15)',
    },

    // ========== Input ==========
    Input: {
      borderRadius: designTokens.borderRadius.lg,
      controlHeight: 36,
      controlHeightLG: 44,
      controlHeightSM: 28,
      paddingBlock: designTokens.spacing.sm,
      paddingInline: designTokens.spacing.md,
    },

    // ========== Select ==========
    Select: {
      borderRadius: designTokens.borderRadius.lg,
      controlHeight: 36,
      controlHeightLG: 44,
      controlHeightSM: 28,
    },

    // ========== Card ==========
    Card: {
      borderRadiusLG: designTokens.borderRadius.xl,
      paddingLG: designTokens.spacing.lg,
      boxShadowTertiary: designTokens.shadows.md,
    },

    // ========== Modal ==========
    Modal: {
      borderRadiusLG: designTokens.borderRadius['2xl'],
      paddingContentHorizontalLG: designTokens.spacing.lg,
      paddingMD: designTokens.spacing.lg,
    },

    // ========== Table ==========
    Table: {
      borderRadius: designTokens.borderRadius.lg,
      headerBg: designTokens.colors.background.secondary,
      headerColor: designTokens.colors.text.primary,
      headerSortActiveBg: designTokens.colors.background.hover,
      rowHoverBg: designTokens.colors.background.hover,
      cellPaddingBlock: designTokens.spacing.md,
      cellPaddingInline: designTokens.spacing.md,
    },

    // ========== Menu ==========
    Menu: {
      itemBorderRadius: designTokens.borderRadius.lg,
      itemHeight: 40,
      itemPaddingInline: designTokens.spacing.md,
      iconSize: 18,
    },

    // ========== Dropdown ==========
    Dropdown: {
      borderRadiusLG: designTokens.borderRadius.lg,
      paddingBlock: designTokens.spacing.sm,
    },

    // ========== Tag ==========
    Tag: {
      borderRadiusSM: designTokens.borderRadius.md,
      defaultBg: designTokens.colors.background.tertiary,
      defaultColor: designTokens.colors.text.secondary,
    },

    // ========== Badge ==========
    Badge: {
      dotSize: 8,
      indicatorHeight: 20,
    },

    // ========== Message ==========
    Message: {
      contentBg: designTokens.colors.white,
      borderRadiusLG: designTokens.borderRadius.xl,
    },

    // ========== Notification ==========
    Notification: {
      width: 384,
      borderRadiusLG: designTokens.borderRadius.xl,
    },

    // ========== Tooltip ==========
    Tooltip: {
      borderRadius: designTokens.borderRadius.md,
      paddingXS: designTokens.spacing.sm,
    },

    // ========== Progress ==========
    Progress: {
      defaultColor: designTokens.colors.primary.default,
      remainingColor: designTokens.colors.background.tertiary,
      circleTextFontSize: `${designTokens.typography.fontSize.lg}px`,
    },

    // ========== Statistic ==========
    Statistic: {
      titleFontSize: designTokens.typography.fontSize.sm,
      contentFontSize: designTokens.typography.fontSize['2xl'],
    },

    // ========== Timeline ==========
    Timeline: {
      dotBorderWidth: 3,
      dotBg: designTokens.colors.white,
      tailColor: designTokens.colors.border.default,
    },

    // ========== Steps ==========
    Steps: {
      dotSize: 32,
      iconSize: 16,
      iconFontSize: 14,
    },

    // ========== Upload ==========
    Upload: {
      actionsColor: designTokens.colors.text.secondary,
    },

    // ========== Tabs ==========
    Tabs: {
      itemActiveColor: designTokens.colors.primary.default,
      itemHoverColor: designTokens.colors.primary.hover,
      itemSelectedColor: designTokens.colors.primary.default,
      inkBarColor: designTokens.colors.primary.default,
      titleFontSize: designTokens.typography.fontSize.md,
    },

    // ========== Breadcrumb ==========
    Breadcrumb: {
      fontSize: designTokens.typography.fontSize.sm,
      iconFontSize: 14,
      linkColor: designTokens.colors.text.secondary,
      linkHoverColor: designTokens.colors.primary.default,
      separatorColor: designTokens.colors.text.tertiary,
    },

    // ========== Pagination ==========
    Pagination: {
      itemSize: 32,
      itemSizeSM: 24,
      borderRadius: designTokens.borderRadius.md,
    },

    // ========== Spin ==========
    Spin: {
      dotSize: 20,
      dotSizeSM: 14,
      dotSizeLG: 32,
    },

    // ========== Alert ==========
    Alert: {
      borderRadiusLG: designTokens.borderRadius.lg,
      withDescriptionPadding: `${designTokens.spacing.md}px ${designTokens.spacing.lg}px`,
    },

    // ========== Divider ==========
    Divider: {
      marginLG: designTokens.spacing.lg,
      textPaddingInline: designTokens.spacing.md,
    },

    // ========== Skeleton ==========
    Skeleton: {
      borderRadiusSM: designTokens.borderRadius.md,
    },

    // ========== Result ==========
    Result: {
      titleFontSize: designTokens.typography.fontSize['2xl'],
      subtitleFontSize: designTokens.typography.fontSize.sm,
      iconFontSize: 72,
    },
  },
}

// Dark theme configuration (for future use)
export const darkTheme: ThemeConfig = {
  ...antdTheme,
  token: {
    ...antdTheme.token,
    colorBgBase: designTokens.colors.gray[900],
    colorTextBase: designTokens.colors.gray[100],
    colorBorder: designTokens.colors.gray[700],
    colorBorderSecondary: designTokens.colors.gray[800],
  },
}

export default antdTheme
