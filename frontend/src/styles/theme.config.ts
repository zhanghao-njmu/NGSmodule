/**
 * Ant Design Theme Configuration
 * Customizes Ant Design components to match our design system
 *
 * @see https://ant.design/docs/react/customize-theme
 */

import type { ThemeConfig } from 'antd'
import { DesignTokens } from './design-tokens'

/**
 * Light theme configuration
 * Default theme for the application
 */
export const lightTheme: ThemeConfig = {
  token: {
    // ===== Colors =====
    colorPrimary: DesignTokens.colors.primary.main,
    colorSuccess: DesignTokens.colors.success.main,
    colorWarning: DesignTokens.colors.warning.main,
    colorError: DesignTokens.colors.error.main,
    colorInfo: DesignTokens.colors.info.main,

    colorTextBase: DesignTokens.colors.text.primary,
    colorBgBase: DesignTokens.colors.background.primary,

    // ===== Typography =====
    fontFamily: DesignTokens.typography.fontFamily.sans,
    fontFamilyCode: DesignTokens.typography.fontFamily.mono,
    fontSize: DesignTokens.typography.fontSize.sm,
    fontSizeHeading1: DesignTokens.typography.fontSize['3xl'],
    fontSizeHeading2: DesignTokens.typography.fontSize['2xl'],
    fontSizeHeading3: DesignTokens.typography.fontSize.xl,
    fontSizeHeading4: DesignTokens.typography.fontSize.lg,
    fontSizeHeading5: DesignTokens.typography.fontSize.md,

    lineHeight: DesignTokens.typography.lineHeight.normal,

    // ===== Border Radius =====
    borderRadius: DesignTokens.radius.md,
    borderRadiusLG: DesignTokens.radius.lg,
    borderRadiusSM: DesignTokens.radius.sm,
    borderRadiusXS: DesignTokens.radius.sm,

    // ===== Spacing =====
    padding: DesignTokens.spacing.md,
    paddingLG: DesignTokens.spacing.lg,
    paddingSM: DesignTokens.spacing.sm,
    paddingXS: DesignTokens.spacing.xs,
    paddingXXS: DesignTokens.spacing.xs / 2,

    margin: DesignTokens.spacing.md,
    marginLG: DesignTokens.spacing.lg,
    marginSM: DesignTokens.spacing.sm,
    marginXS: DesignTokens.spacing.xs,
    marginXXS: DesignTokens.spacing.xs / 2,

    // ===== Shadows =====
    boxShadow: DesignTokens.shadow.md,
    boxShadowSecondary: DesignTokens.shadow.sm,

    // ===== Motion =====
    motionDurationFast: `${DesignTokens.duration.fast}ms`,
    motionDurationMid: `${DesignTokens.duration.base}ms`,
    motionDurationSlow: `${DesignTokens.duration.slow}ms`,

    motionEaseInOut: 'cubic-bezier(0.4, 0, 0.2, 1)',
    motionEaseOut: 'cubic-bezier(0, 0, 0.2, 1)',
    motionEaseIn: 'cubic-bezier(0.4, 0, 1, 1)',

    // ===== Z-Index =====
    zIndexPopupBase: DesignTokens.zIndex.dropdown,
    zIndexBase: DesignTokens.zIndex.dropdown,
  },

  components: {
    // ===== Button =====
    Button: {
      borderRadius: DesignTokens.radius.md,
      borderRadiusLG: DesignTokens.radius.lg,
      borderRadiusSM: DesignTokens.radius.sm,
      fontWeight: DesignTokens.typography.fontWeight.medium,
      primaryShadow: '0 2px 4px rgba(37, 99, 235, 0.2)',
      defaultShadow: '0 1px 2px rgba(0, 0, 0, 0.05)',
    },

    // ===== Card =====
    Card: {
      borderRadiusLG: DesignTokens.radius.lg,
      boxShadowTertiary: DesignTokens.shadow.md,
      paddingLG: DesignTokens.spacing.lg,
      headerHeight: 48,
      headerFontSize: DesignTokens.typography.fontSize.lg,
      headerFontSizeSM: DesignTokens.typography.fontSize.md,
    },

    // ===== Table =====
    Table: {
      borderRadius: DesignTokens.radius.lg,
      headerBg: DesignTokens.colors.gray[50],
      headerColor: DesignTokens.colors.text.primary,
      headerSortActiveBg: DesignTokens.colors.gray[100],
      bodySortBg: DesignTokens.colors.primary.lighter,
      rowHoverBg: DesignTokens.colors.primary.lighter,
      headerBorderRadius: DesignTokens.radius.lg,
    },

    // ===== Modal =====
    Modal: {
      borderRadiusLG: DesignTokens.radius.lg,
      headerBg: DesignTokens.colors.background.primary,
      contentBg: DesignTokens.colors.background.primary,
      footerBg: DesignTokens.colors.background.primary,
      bodyPadding: DesignTokens.spacing.lg,
      footerPadding: `${DesignTokens.spacing.md}px ${DesignTokens.spacing.lg}px`,
    },

    // ===== Input =====
    Input: {
      borderRadius: DesignTokens.radius.md,
      borderRadiusLG: DesignTokens.radius.lg,
      borderRadiusSM: DesignTokens.radius.sm,
      hoverBorderColor: DesignTokens.colors.primary.main,
      activeBorderColor: DesignTokens.colors.primary.main,
      activeShadow: `0 0 0 2px ${DesignTokens.colors.primary.light}`,
    },

    // ===== Select =====
    Select: {
      borderRadius: DesignTokens.radius.md,
      borderRadiusLG: DesignTokens.radius.lg,
      borderRadiusSM: DesignTokens.radius.sm,
    },

    // ===== Tag =====
    Tag: {
      borderRadiusSM: DesignTokens.radius.md,
      fontSizeSM: DesignTokens.typography.fontSize.xs,
    },

    // ===== Progress =====
    Progress: {
      lineBorderRadius: DesignTokens.radius.full,
      circleTextFontSize: `${DesignTokens.typography.fontSize.xl}px`,
    },

    // ===== Message =====
    Message: {
      contentBg: DesignTokens.colors.background.primary,
      borderRadius: DesignTokens.radius.lg,
    },

    // ===== Notification =====
    Notification: {
      borderRadiusLG: DesignTokens.radius.lg,
    },

    // ===== Menu =====
    Menu: {
      itemBorderRadius: DesignTokens.radius.md,
      itemSelectedBg: DesignTokens.colors.primary.light,
      itemSelectedColor: DesignTokens.colors.primary.main,
      itemActiveBg: DesignTokens.colors.primary.lighter,
      subMenuItemBorderRadius: DesignTokens.radius.md,
    },

    // ===== Tabs =====
    Tabs: {
      horizontalItemPadding: `${DesignTokens.spacing.md}px ${DesignTokens.spacing.lg}px`,
      titleFontSize: DesignTokens.typography.fontSize.md,
      inkBarColor: DesignTokens.colors.primary.main,
      itemActiveColor: DesignTokens.colors.primary.main,
      itemHoverColor: DesignTokens.colors.primary.hover,
      itemSelectedColor: DesignTokens.colors.primary.main,
    },

    // ===== Badge =====
    Badge: {
      dotSize: 6,
      textFontSize: DesignTokens.typography.fontSize.xs,
      textFontSizeSM: DesignTokens.typography.fontSize.xs,
    },

    // ===== Dropdown =====
    Dropdown: {
      borderRadiusLG: DesignTokens.radius.lg,
      borderRadiusSM: DesignTokens.radius.md,
    },

    // ===== Tooltip =====
    Tooltip: {
      borderRadius: DesignTokens.radius.md,
    },

    // ===== Popover =====
    Popover: {
      borderRadiusLG: DesignTokens.radius.lg,
    },

    // ===== Alert =====
    Alert: {
      borderRadiusLG: DesignTokens.radius.lg,
      withDescriptionPadding: `${DesignTokens.spacing.md}px ${DesignTokens.spacing.lg}px`,
    },

    // ===== Statistic =====
    Statistic: {
      titleFontSize: DesignTokens.typography.fontSize.sm,
      contentFontSize: DesignTokens.typography.fontSize['2xl'],
    },

    // ===== Timeline =====
    Timeline: {
      tailColor: DesignTokens.colors.border.default,
      dotBorderWidth: 2,
    },

    // ===== Steps =====
    Steps: {
      dotSize: 32,
      iconSize: 32,
      iconSizeSM: 24,
    },

    // ===== Skeleton =====
    Skeleton: {
      borderRadiusSM: DesignTokens.radius.md,
    },

    // ===== Spin =====
    Spin: {
      dotSize: 20,
      dotSizeSM: 14,
      dotSizeLG: 32,
    },

    // ===== Upload =====
    Upload: {
      actionsColor: DesignTokens.colors.text.tertiary,
    },

    // ===== Drawer =====
    Drawer: {
      footerPaddingBlock: DesignTokens.spacing.md,
      footerPaddingInline: DesignTokens.spacing.lg,
    },

    // ===== Switch =====
    Switch: {
      trackMinWidth: 44,
      trackHeight: 22,
      trackPadding: 2,
      innerMinMargin: 4,
      innerMaxMargin: 24,
    },

    // ===== DatePicker =====
    DatePicker: {
      borderRadius: DesignTokens.radius.md,
      borderRadiusLG: DesignTokens.radius.lg,
      borderRadiusSM: DesignTokens.radius.sm,
    },

    // ===== Pagination =====
    Pagination: {
      borderRadius: DesignTokens.radius.md,
      itemActiveBg: DesignTokens.colors.primary.main,
      itemLinkBg: DesignTokens.colors.background.primary,
      itemBg: DesignTokens.colors.background.primary,
    },
  },

  // ===== Algorithm =====
  // Uses Ant Design's default light algorithm
  algorithm: undefined,
}

/**
 * Dark theme configuration
 * Theme for dark mode
 */
export const darkTheme: ThemeConfig = {
  ...lightTheme,

  token: {
    ...lightTheme.token,

    // Override colors for dark mode
    colorTextBase: DesignTokens.colors.gray[100],
    colorBgBase: DesignTokens.colors.gray[900],

    colorBgContainer: DesignTokens.colors.gray[800],
    colorBgElevated: DesignTokens.colors.gray[800],
    colorBgLayout: DesignTokens.colors.gray[900],

    colorBorder: DesignTokens.colors.gray[700],
    colorBorderSecondary: DesignTokens.colors.gray[700],
  },

  components: {
    ...lightTheme.components,

    Card: {
      ...lightTheme.components?.Card,
      colorBgContainer: DesignTokens.colors.gray[800],
    },

    Table: {
      ...lightTheme.components?.Table,
      headerBg: DesignTokens.colors.gray[800],
      bodySortBg: DesignTokens.colors.gray[700],
      rowHoverBg: DesignTokens.colors.gray[700],
    },

    Menu: {
      ...lightTheme.components?.Menu,
      itemBg: 'transparent',
      itemSelectedBg: DesignTokens.colors.gray[700],
      itemActiveBg: DesignTokens.colors.gray[800],
      subMenuItemBg: 'transparent',
    },
  },
}

/**
 * Get theme configuration based on mode
 */
export const getTheme = (mode: 'light' | 'dark' = 'light'): ThemeConfig => {
  return mode === 'dark' ? darkTheme : lightTheme
}

export default lightTheme
