/**
 * Ant Design Theme Configuration
 * NGSmodule - Modern Bioinformatics Platform
 */

import type { ThemeConfig } from 'antd'

export const theme: ThemeConfig = {
  token: {
    // ========== Colors ==========
    colorPrimary: '#2563eb', // Scientific Blue
    colorSuccess: '#10b981', // Green
    colorWarning: '#f59e0b', // Amber
    colorError: '#ef4444', // Red
    colorInfo: '#06b6d4', // Cyan

    // ========== Background ==========
    colorBgBase: '#ffffff',
    colorBgContainer: '#ffffff',
    colorBgElevated: '#ffffff',
    colorBgLayout: '#f9fafb',

    // ========== Text ==========
    colorText: '#111827',
    colorTextSecondary: '#6b7280',
    colorTextTertiary: '#9ca3af',
    colorTextQuaternary: '#d1d5db',

    // ========== Border ==========
    colorBorder: '#e5e7eb',
    colorBorderSecondary: '#f3f4f6',

    // ========== Typography ==========
    fontSize: 14,
    fontSizeHeading1: 30,
    fontSizeHeading2: 24,
    fontSizeHeading3: 20,
    fontSizeHeading4: 18,
    fontSizeHeading5: 16,
    fontWeightStrong: 600,

    // ========== Spacing ==========
    sizeStep: 4,
    sizeUnit: 4,
    padding: 16,
    paddingLG: 24,
    paddingXL: 32,
    margin: 16,
    marginLG: 24,
    marginXL: 32,

    // ========== Border Radius ==========
    borderRadius: 8,
    borderRadiusLG: 12,
    borderRadiusSM: 4,

    // ========== Line Height ==========
    lineHeight: 1.5,
    lineHeightHeading1: 1.25,
    lineHeightHeading2: 1.25,

    // ========== Control ==========
    controlHeight: 36,
    controlHeightLG: 40,
    controlHeightSM: 28,

    // ========== Box Shadow ==========
    boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
    boxShadowSecondary: '0 1px 2px 0 rgba(0, 0, 0, 0.05)',

    // ========== Motion ==========
    motionDurationFast: '0.15s',
    motionDurationMid: '0.25s',
    motionDurationSlow: '0.35s',
  },

  components: {
    // ========== Button ==========
    Button: {
      primaryShadow: '0 2px 4px rgba(37, 99, 235, 0.2)',
      controlHeight: 36,
      controlHeightLG: 40,
      controlHeightSM: 28,
      fontWeight: 500,
      borderRadius: 8,
      paddingContentHorizontal: 16,
    },

    // ========== Card ==========
    Card: {
      borderRadiusLG: 12,
      boxShadow: '0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)',
      headerFontSize: 16,
      headerFontSizeSM: 14,
      headerHeight: 56,
      headerHeightSM: 48,
      paddingLG: 24,
      padding: 20,
    },

    // ========== Table ==========
    Table: {
      headerBg: '#f9fafb',
      headerColor: '#111827',
      headerSplitColor: '#e5e7eb',
      rowHoverBg: '#eff6ff',
      borderRadius: 8,
      cellPaddingBlock: 12,
      cellPaddingInline: 16,
      headerBorderRadius: 8,
    },

    // ========== Input ==========
    Input: {
      controlHeight: 36,
      controlHeightLG: 40,
      controlHeightSM: 28,
      borderRadius: 8,
      paddingBlock: 8,
      paddingInline: 12,
      activeShadow: '0 0 0 2px #dbeafe',
    },

    // ========== Select ==========
    Select: {
      controlHeight: 36,
      controlHeightLG: 40,
      controlHeightSM: 28,
      borderRadius: 8,
      optionSelectedBg: '#eff6ff',
      optionActiveBg: '#f9fafb',
    },

    // ========== Modal ==========
    Modal: {
      borderRadiusLG: 12,
      headerBg: '#ffffff',
      contentBg: '#ffffff',
      footerBg: '#ffffff',
      titleFontSize: 18,
      titleLineHeight: 1.5,
    },

    // ========== Tag ==========
    Tag: {
      borderRadiusSM: 6,
      fontSizeSM: 12,
    },

    // ========== Menu ==========
    Menu: {
      itemBorderRadius: 8,
      itemHeight: 40,
      itemMarginInline: 8,
      itemPaddingInline: 16,
      itemSelectedBg: '#eff6ff',
      itemSelectedColor: '#2563eb',
      itemHoverBg: '#f9fafb',
    },

    // ========== Tabs ==========
    Tabs: {
      itemActiveColor: '#2563eb',
      itemHoverColor: '#1d4ed8',
      itemSelectedColor: '#2563eb',
      cardBg: '#ffffff',
      cardHeight: 44,
      titleFontSize: 14,
      titleFontSizeLG: 16,
      inkBarColor: '#2563eb',
    },

    // ========== Progress ==========
    Progress: {
      remainingColor: '#e5e7eb',
      defaultColor: '#2563eb',
    },

    // ========== Badge ==========
    Badge: {
      indicatorHeight: 20,
      indicatorHeightSM: 16,
      dotSize: 6,
      textFontSize: 12,
      textFontSizeSM: 10,
    },

    // ========== Alert ==========
    Alert: {
      borderRadiusLG: 8,
      defaultPadding: '12px 16px',
      withDescriptionPadding: '16px 20px',
      withDescriptionIconSize: 24,
    },

    // ========== Message ==========
    Message: {
      contentBg: '#ffffff',
      contentPadding: '12px 16px',
      borderRadiusLG: 8,
    },

    // ========== Notification ==========
    Notification: {
      width: 384,
      borderRadiusLG: 12,
      padding: 20,
      paddingLG: 24,
    },

    // ========== Dropdown ==========
    Dropdown: {
      borderRadiusLG: 8,
      paddingBlock: 8,
      controlPaddingHorizontal: 12,
    },

    // ========== Tooltip ==========
    Tooltip: {
      borderRadius: 6,
      fontSize: 12,
    },

    // ========== Statistic ==========
    Statistic: {
      titleFontSize: 14,
      contentFontSize: 24,
    },

    // ========== Pagination ==========
    Pagination: {
      itemSize: 32,
      itemSizeSM: 24,
      borderRadius: 6,
      itemActiveBg: '#2563eb',
      itemLinkBg: '#ffffff',
      itemBg: '#ffffff',
    },

    // ========== Switch ==========
    Switch: {
      trackHeight: 22,
      trackHeightSM: 16,
      trackMinWidth: 44,
      trackMinWidthSM: 28,
      handleSize: 18,
      handleSizeSM: 12,
      innerMinMargin: 4,
      innerMaxMargin: 24,
    },

    // ========== Divider ==========
    Divider: {
      marginLG: 24,
      margin: 16,
      orientationMargin: 0.05,
      verticalMarginInline: 8,
      colorSplit: '#e5e7eb',
    },

    // ========== Steps ==========
    Steps: {
      navArrowColor: '#2563eb',
      dotCurrentSize: 10,
      dotSize: 8,
      titleLineHeight: 32,
      iconSize: 32,
      iconSizeSM: 24,
    },

    // ========== Upload ==========
    Upload: {
      actionsColor: '#2563eb',
      borderRadiusLG: 8,
      paddingLG: 24,
    },
  },

  // ========== Algorithm ==========
  algorithm: undefined, // 使用默认算法，可以切换为 theme.darkAlgorithm
}

export default theme
