# NGSmodule 设计系统

## 概述

NGSmodule 采用现代化、科学化的设计系统，专为生物信息学研究人员设计，强调可读性、专业性和易用性。

## 设计原则

### 1. 科学专业 (Scientific & Professional)
- 使用科学蓝作为主色调，传达专业和可信赖感
- 清晰的数据展示，适合复杂的生物信息学数据
- 专业的配色方案，避免过度花哨

### 2. 清晰简洁 (Clear & Concise)
- 简洁的布局，突出核心功能
- 明确的层级关系
- 充足的留白空间

### 3. 高效易用 (Efficient & Usable)
- 快速的操作流程
- 直观的交互反馈
- 减少认知负担

### 4. 一致性 (Consistency)
- 统一的视觉语言
- 一致的交互模式
- 可预测的行为

---

## 颜色系统

### 主色 (Primary Color)
**科学蓝** - 专业、可信赖

```css
--color-primary: #2563eb
--color-primary-hover: #1d4ed8
--color-primary-active: #1e40af
--color-primary-light: #dbeafe
--color-primary-lighter: #eff6ff
```

**用途**:
- 主要按钮
- 链接
- 重要信息高亮
- 进行中的任务状态

### 功能色 (Functional Colors)

#### 成功 (Success) - 绿色
```css
--color-success: #10b981
```
- 完成状态
- 成功通知
- 正向操作反馈

#### 警告 (Warning) - 琥珀色
```css
--color-warning: #f59e0b
```
- 警告信息
- 暂停状态
- 需要注意的操作

#### 错误 (Error) - 红色
```css
--color-error: #ef4444
```
- 失败状态
- 错误信息
- 危险操作

#### 信息 (Info) - 青色
```css
--color-info: #06b6d4
```
- 提示信息
- 辅助说明

### 中性色 (Neutral Colors)

```css
--color-gray-50: #f9fafb   /* 背景色 */
--color-gray-100: #f3f4f6  /* 次级背景 */
--color-gray-200: #e5e7eb  /* 边框 */
--color-gray-300: #d1d5db  /* 分隔线 */
--color-gray-400: #9ca3af  /* 禁用文本 */
--color-gray-500: #6b7280  /* 辅助文本 */
--color-gray-600: #4b5563  /* 次要文本 */
--color-gray-700: #374151  /* 标题 */
--color-gray-800: #1f2937  /* 强调文本 */
--color-gray-900: #111827  /* 主要文本 */
```

### 状态色映射

```typescript
{
  pending: gray-400,      // 等待中
  running: primary,       // 运行中
  completed: success,     // 已完成
  failed: error,          // 失败
  paused: warning,        // 暂停
  cancelled: gray-500,    // 已取消
}
```

### 文件类型色映射

```typescript
{
  fastq: primary (#2563eb),    // FASTQ 文件
  bam: success (#10b981),       // BAM 文件
  vcf: purple (#722ed1),        // VCF 文件
  sam: warning (#f59e0b),       // SAM 文件
  bed: info (#06b6d4),          // BED 文件
  gff: magenta (#eb2f96),       // GFF 文件
}
```

---

## 排版系统

### 字体家族

```css
--font-family-sans: -apple-system, BlinkMacSystemFont, 'Segoe UI',
                    Roboto, 'Helvetica Neue', Arial, sans-serif
--font-family-mono: 'Fira Code', 'Courier New', monospace
```

### 字号标尺 (Type Scale)

| 名称 | 大小 | 用途 |
|------|------|------|
| xs | 12px | 辅助说明、时间戳 |
| sm | 14px | 正文、表单标签 |
| md | 16px | 强调文本、按钮 |
| lg | 18px | 小标题 |
| xl | 20px | 卡片标题 |
| 2xl | 24px | 页面标题 |
| 3xl | 30px | 大标题 |
| 4xl | 36px | 特大标题、数字展示 |

### 字重 (Font Weight)

| 名称 | 值 | 用途 |
|------|-----|------|
| normal | 400 | 正文 |
| medium | 500 | 按钮、标签 |
| semibold | 600 | 小标题 |
| bold | 700 | 强调标题 |

### 行高 (Line Height)

```css
--line-height-tight: 1.25      /* 标题 */
--line-height-normal: 1.5      /* 正文 */
--line-height-relaxed: 1.75    /* 长文本 */
```

---

## 间距系统

基于 **4px 网格系统**

```css
--spacing-xs: 4px     /* 0.25rem */
--spacing-sm: 8px     /* 0.5rem */
--spacing-md: 16px    /* 1rem */
--spacing-lg: 24px    /* 1.5rem */
--spacing-xl: 32px    /* 2rem */
--spacing-2xl: 48px   /* 3rem */
--spacing-3xl: 64px   /* 4rem */
```

### 使用指南

- **xs (4px)**: 图标间距、紧密元素
- **sm (8px)**: 表单元素内边距、列表项间距
- **md (16px)**: 卡片内边距、段落间距
- **lg (24px)**: 模块间距、卡片外边距
- **xl (32px)**: 页面主要区块间距
- **2xl (48px)**: 大区块间距
- **3xl (64px)**: 页面顶部/底部留白

---

## 圆角系统

```css
--radius-sm: 2px      /* 小元素: Badge */
--radius-md: 4px      /* 表单元素内部 */
--radius-lg: 8px      /* 按钮、输入框、卡片 */
--radius-xl: 12px     /* 大卡片、Modal */
--radius-2xl: 16px    /* 特殊容器 */
--radius-full: 9999px /* 圆形: Avatar, Tag */
```

---

## 阴影系统

```css
--shadow-sm: 0 1px 2px rgba(0, 0, 0, 0.05)
  /* 微弱阴影: 悬浮按钮 */

--shadow-md: 0 4px 6px -1px rgba(0, 0, 0, 0.1)
  /* 标准阴影: 卡片、下拉菜单 */

--shadow-lg: 0 10px 15px -3px rgba(0, 0, 0, 0.1)
  /* 明显阴影: 弹出层 */

--shadow-xl: 0 20px 25px -5px rgba(0, 0, 0, 0.1)
  /* 强阴影: Modal */

--shadow-2xl: 0 25px 50px -12px rgba(0, 0, 0, 0.25)
  /* 极强阴影: 特殊浮动元素 */
```

---

## 动效系统

### 过渡时长

```css
--duration-fast: 150ms     /* 快速: Hover, Focus */
--duration-base: 250ms     /* 标准: 状态切换 */
--duration-slow: 350ms     /* 缓慢: Modal打开/关闭 */
```

### 缓动函数

```css
cubic-bezier(0.4, 0, 0.2, 1)  /* 标准缓动 */
```

### 应用场景

- **Hover**: fast (150ms)
- **Focus**: fast (150ms)
- **Toggle**: base (250ms)
- **Collapse/Expand**: slow (350ms)
- **Modal**: slow (350ms)
- **Page Transition**: base (250ms)

---

## 组件规范

### 按钮 (Button)

**尺寸**:
- Small: 28px
- Default: 36px
- Large: 44px

**边距**:
- Horizontal Padding: 16px
- Icon Spacing: 8px

**圆角**: 8px

### 输入框 (Input)

**尺寸**:
- Small: 28px
- Default: 36px
- Large: 44px

**边距**:
- Horizontal Padding: 12px
- Vertical Padding: 8px

**圆角**: 8px

### 卡片 (Card)

**内边距**: 24px
**圆角**: 12px
**阴影**: shadow-md

### 表格 (Table)

**表头背景**: gray-50
**行高亮**: primary-lighter
**单元格内边距**: 12px (vertical) × 16px (horizontal)
**圆角**: 8px

---

## 布局系统

### 断点 (Breakpoints)

```css
--breakpoint-xs: 480px    /* 手机 */
--breakpoint-sm: 576px    /* 大手机 */
--breakpoint-md: 768px    /* 平板 */
--breakpoint-lg: 992px    /* 桌面 */
--breakpoint-xl: 1200px   /* 大桌面 */
--breakpoint-2xl: 1600px  /* 超大桌面 */
```

### 布局尺寸

```css
--sidebar-width: 240px
--sidebar-collapsed-width: 64px
--header-height: 64px
--footer-height: 48px
```

### 容器最大宽度

- 内容区域: 1200px
- 宽屏内容: 1600px
- 全宽: 100%

---

## 图标规范

### 尺寸

- Small: 14px
- Default: 16px
- Medium: 18px
- Large: 20px
- Extra Large: 24px

### 使用原则

- 与文字对齐时，图标大小应匹配字号
- 操作按钮图标: 16px - 18px
- 卡片/列表图标: 18px - 20px
- 大型展示图标: 48px - 72px

---

## 无障碍 (Accessibility)

### 颜色对比度

- 正文文本: 至少 4.5:1
- 大文本(18px+): 至少 3:1
- 交互元素: 至少 3:1

### 焦点指示

- 明确的 Focus 状态
- 使用 outline 或 box-shadow
- 颜色: primary-light (#dbeafe)

### 键盘导航

- 所有交互元素可键盘访问
- Tab 顺序符合逻辑
- 提供快捷键提示

---

## TypeScript 使用

### 导入设计令牌

```typescript
import { designTokens } from '@/config/design-tokens'

// 使用颜色
const primaryColor = designTokens.colors.primary.default

// 使用间距
const padding = designTokens.spacing.md

// 使用状态色
import { statusColors } from '@/config/design-tokens'
const color = statusColors.completed  // #10b981
```

### 类型支持

```typescript
import type { DesignTokens, ColorKey, SpacingKey } from '@/config/design-tokens'
```

---

## 最佳实践

### ✅ 推荐

- 优先使用 CSS 变量而非硬编码值
- 保持组件样式一致性
- 遵循 4px 网格系统
- 使用语义化的颜色名称
- 为交互元素提供明确反馈

### ❌ 避免

- 随意使用非系统颜色
- 破坏间距规律
- 过度使用动画
- 忽略暗色模式兼容性
- 忽略无障碍需求

---

## 暗色模式 (Dark Mode)

当前已预留暗色模式变量，可通过 `[data-theme='dark']` 切换：

```css
[data-theme='dark'] {
  --bg-primary: #111827;
  --text-primary: #f3f4f6;
  --border-color: #374151;
  /* ... */
}
```

---

## 参考资源

- Ant Design 5: https://ant.design/
- Tailwind CSS: https://tailwindcss.com/ (颜色参考)
- Material Design: https://material.io/ (设计原则参考)

---

**版本**: 1.0.0
**最后更新**: 2025-01-22
