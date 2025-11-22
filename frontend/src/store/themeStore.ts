import { create } from 'zustand'
import { persist } from 'zustand/middleware'

type ThemeMode = 'light' | 'dark'

interface ThemeState {
  mode: ThemeMode
  setMode: (mode: ThemeMode) => void
  toggleMode: () => void
}

/**
 * Theme Store
 * Manages the application theme (light/dark mode)
 * Persists theme preference to localStorage
 */
export const useThemeStore = create<ThemeState>()(
  persist(
    (set, get) => ({
      mode: 'light',

      setMode: (mode: ThemeMode) => {
        set({ mode })
        // Update CSS data attribute for dark mode styles
        document.documentElement.setAttribute('data-theme', mode)
      },

      toggleMode: () => {
        const currentMode = get().mode
        const newMode = currentMode === 'light' ? 'dark' : 'light'
        get().setMode(newMode)
      },
    }),
    {
      name: 'ngsmodule-theme',
      onRehydrateStorage: () => (state) => {
        // Apply theme on initial load
        if (state) {
          document.documentElement.setAttribute('data-theme', state.mode)
        }
      },
    },
  ),
)

/**
 * Custom hook for theme management
 * Provides theme mode and toggle functionality
 *
 * @example
 * const { mode, toggleMode } = useTheme()
 */
export const useTheme = () => {
  const { mode, setMode, toggleMode } = useThemeStore()

  return {
    mode,
    isDark: mode === 'dark',
    isLight: mode === 'light',
    setMode,
    toggleMode,
  }
}
