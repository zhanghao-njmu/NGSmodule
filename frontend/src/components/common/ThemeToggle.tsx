import { Switch } from 'antd'
import { BulbOutlined, BulbFilled } from '@ant-design/icons'
import { useTheme } from '@/store/themeStore'
import './ThemeToggle.css'

interface ThemeToggleProps {
  /**
   * Display mode: icon-only or full (with text)
   */
  mode?: 'icon' | 'full'
  /**
   * Size of the toggle
   */
  size?: 'small' | 'default'
}

/**
 * Theme Toggle Component
 * Allows users to switch between light and dark mode
 *
 * @example
 * // Icon only
 * <ThemeToggle mode="icon" />
 *
 * // With text
 * <ThemeToggle mode="full" />
 */
export const ThemeToggle: React.FC<ThemeToggleProps> = ({ mode = 'icon', size = 'default' }) => {
  const { isDark, toggleMode } = useTheme()

  if (mode === 'icon') {
    return (
      <div className="theme-toggle-icon" onClick={toggleMode}>
        {isDark ? (
          <BulbFilled className="theme-icon theme-icon-dark" style={{ fontSize: size === 'small' ? 18 : 20 }} />
        ) : (
          <BulbOutlined className="theme-icon theme-icon-light" style={{ fontSize: size === 'small' ? 18 : 20 }} />
        )}
      </div>
    )
  }

  return (
    <div className="theme-toggle-full">
      <BulbOutlined className="theme-icon-inline" />
      <span className="theme-label">深色模式</span>
      <Switch checked={isDark} onChange={toggleMode} size={size} checkedChildren="开" unCheckedChildren="关" />
    </div>
  )
}
