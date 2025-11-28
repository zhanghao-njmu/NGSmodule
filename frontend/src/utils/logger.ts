/**
 * Logger utility for frontend
 *
 * Provides consistent logging with environment-aware behavior:
 * - In development: logs to console
 * - In production: logs only warnings and errors
 *
 * Usage:
 *   import { logger } from '@/utils/logger'
 *   logger.info('message')
 *   logger.warn('warning')
 *   logger.error('error', errorObject)
 */

const isDevelopment = import.meta.env.DEV

type LogLevel = 'debug' | 'info' | 'warn' | 'error'

interface Logger {
  debug: (...args: unknown[]) => void
  info: (...args: unknown[]) => void
  warn: (...args: unknown[]) => void
  error: (...args: unknown[]) => void
  log: (...args: unknown[]) => void
}

function createLogger(): Logger {
  const shouldLog = (level: LogLevel): boolean => {
    if (isDevelopment) {
      return true
    }
    // In production, only log warnings and errors
    return level === 'warn' || level === 'error'
  }

  return {
    debug: (...args: unknown[]) => {
      if (shouldLog('debug')) {
        console.debug('[DEBUG]', ...args)
      }
    },
    info: (...args: unknown[]) => {
      if (shouldLog('info')) {
        console.info('[INFO]', ...args)
      }
    },
    warn: (...args: unknown[]) => {
      if (shouldLog('warn')) {
        console.warn('[WARN]', ...args)
      }
    },
    error: (...args: unknown[]) => {
      if (shouldLog('error')) {
        console.error('[ERROR]', ...args)
      }
    },
    // Alias for info
    log: (...args: unknown[]) => {
      if (shouldLog('info')) {
        console.log(...args)
      }
    },
  }
}

export const logger = createLogger()
export default logger
