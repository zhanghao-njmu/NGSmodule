/**
 * Formatting Utilities
 * Centralized formatting functions for consistent display across the application
 */

import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'
import duration from 'dayjs/plugin/duration'

// Extend dayjs with plugins
dayjs.extend(relativeTime)
dayjs.extend(duration)

/**
 * Format file size in bytes to human-readable string
 *
 * @param bytes - File size in bytes
 * @param decimals - Number of decimal places (default: 2)
 * @returns Formatted file size string (e.g., "1.50 MB")
 *
 * @example
 * formatFileSize(1536) // "1.50 KB"
 * formatFileSize(1048576) // "1.00 MB"
 * formatFileSize(0) // "0 B"
 */
export const formatFileSize = (bytes: number, decimals: number = 2): string => {
  if (bytes === 0) return '0 B'
  if (bytes < 0) return 'Invalid size'

  const k = 1024
  const dm = decimals < 0 ? 0 : decimals
  const sizes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

  const i = Math.floor(Math.log(bytes) / Math.log(k))
  const size = bytes / Math.pow(k, i)

  return `${size.toFixed(dm)} ${sizes[i]}`
}

/**
 * Format date to standard datetime format (YYYY-MM-DD HH:mm)
 *
 * @param date - Date string or Date object
 * @param format - Custom format string (default: 'YYYY-MM-DD HH:mm')
 * @returns Formatted date string
 *
 * @example
 * formatDateTime('2024-01-15T10:30:00Z') // "2024-01-15 10:30"
 * formatDateTime(new Date(), 'YYYY/MM/DD') // "2024/01/15"
 */
export const formatDateTime = (
  date: string | Date | null | undefined,
  format: string = 'YYYY-MM-DD HH:mm'
): string => {
  if (!date) return '-'
  return dayjs(date).format(format)
}

/**
 * Format date to relative time (e.g., "2 hours ago", "in 3 days")
 *
 * @param date - Date string or Date object
 * @returns Relative time string
 *
 * @example
 * formatDateRelative('2024-01-15T10:00:00Z') // "2 hours ago"
 * formatDateRelative(new Date()) // "a few seconds ago"
 */
export const formatDateRelative = (date: string | Date | null | undefined): string => {
  if (!date) return '-'
  return dayjs(date).fromNow()
}

/**
 * Format date to long format with seconds (YYYY-MM-DD HH:mm:ss)
 *
 * @param date - Date string or Date object
 * @returns Formatted date string with seconds
 *
 * @example
 * formatDateLong('2024-01-15T10:30:45Z') // "2024-01-15 10:30:45"
 */
export const formatDateLong = (date: string | Date | null | undefined): string => {
  return formatDateTime(date, 'YYYY-MM-DD HH:mm:ss')
}

/**
 * Format date to date only (YYYY-MM-DD)
 *
 * @param date - Date string or Date object
 * @returns Formatted date string without time
 *
 * @example
 * formatDate('2024-01-15T10:30:00Z') // "2024-01-15"
 */
export const formatDate = (date: string | Date | null | undefined): string => {
  return formatDateTime(date, 'YYYY-MM-DD')
}

/**
 * Format time only (HH:mm:ss)
 *
 * @param date - Date string or Date object
 * @returns Formatted time string
 *
 * @example
 * formatTime('2024-01-15T10:30:45Z') // "10:30:45"
 */
export const formatTime = (date: string | Date | null | undefined): string => {
  return formatDateTime(date, 'HH:mm:ss')
}

/**
 * Format number with specified decimal places
 *
 * @param value - Number to format
 * @param decimals - Number of decimal places (default: 2)
 * @param fallback - Fallback value if number is invalid (default: '-')
 * @returns Formatted number string
 *
 * @example
 * formatNumber(1234.5678) // "1234.57"
 * formatNumber(1234.5678, 0) // "1235"
 * formatNumber(null) // "-"
 */
export const formatNumber = (
  value: number | null | undefined,
  decimals: number = 2,
  fallback: string = '-'
): string => {
  if (value === null || value === undefined || isNaN(value)) return fallback
  return value.toFixed(decimals)
}

/**
 * Format number with thousands separator
 *
 * @param value - Number to format
 * @param decimals - Number of decimal places (default: 0)
 * @returns Formatted number with commas
 *
 * @example
 * formatNumberWithCommas(1234567) // "1,234,567"
 * formatNumberWithCommas(1234.5678, 2) // "1,234.57"
 */
export const formatNumberWithCommas = (
  value: number | null | undefined,
  decimals: number = 0
): string => {
  if (value === null || value === undefined || isNaN(value)) return '-'

  const parts = value.toFixed(decimals).split('.')
  parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ',')

  return parts.join('.')
}

/**
 * Format percentage (0-1 to percentage string)
 *
 * @param value - Decimal value (0-1)
 * @param decimals - Number of decimal places (default: 1)
 * @param fallback - Fallback value if invalid (default: '-')
 * @returns Formatted percentage string with % symbol
 *
 * @example
 * formatPercent(0.856) // "85.6%"
 * formatPercent(0.5, 0) // "50%"
 * formatPercent(1) // "100.0%"
 */
export const formatPercent = (
  value: number | null | undefined,
  decimals: number = 1,
  fallback: string = '-'
): string => {
  if (value === null || value === undefined || isNaN(value)) return fallback
  return `${(value * 100).toFixed(decimals)}%`
}

/**
 * Format percentage from ratio (numerator/denominator)
 *
 * @param numerator - Numerator value
 * @param denominator - Denominator value
 * @param decimals - Number of decimal places (default: 1)
 * @returns Formatted percentage string
 *
 * @example
 * formatPercentFromRatio(45, 100) // "45.0%"
 * formatPercentFromRatio(3, 7, 2) // "42.86%"
 */
export const formatPercentFromRatio = (
  numerator: number,
  denominator: number,
  decimals: number = 1
): string => {
  if (denominator === 0) return '-'
  return formatPercent(numerator / denominator, decimals)
}

/**
 * Format duration in seconds to human-readable string
 *
 * @param seconds - Duration in seconds
 * @param format - Format style: 'short' (1h 30m 15s) or 'long' (1 hour 30 minutes 15 seconds)
 * @returns Formatted duration string
 *
 * @example
 * formatDuration(90) // "1m 30s"
 * formatDuration(3665) // "1h 1m 5s"
 * formatDuration(45) // "45s"
 */
export const formatDuration = (seconds: number, format: 'short' | 'long' = 'short'): string => {
  if (seconds < 0) return '-'
  if (seconds === 0) return '0s'

  const hours = Math.floor(seconds / 3600)
  const minutes = Math.floor((seconds % 3600) / 60)
  const secs = seconds % 60

  const parts: string[] = []

  if (format === 'short') {
    if (hours > 0) parts.push(`${hours}h`)
    if (minutes > 0) parts.push(`${minutes}m`)
    if (secs > 0 || parts.length === 0) parts.push(`${secs}s`)
  } else {
    if (hours > 0) parts.push(`${hours} ${hours === 1 ? 'hour' : 'hours'}`)
    if (minutes > 0) parts.push(`${minutes} ${minutes === 1 ? 'minute' : 'minutes'}`)
    if (secs > 0 || parts.length === 0) parts.push(`${secs} ${secs === 1 ? 'second' : 'seconds'}`)
  }

  return parts.join(' ')
}

/**
 * Format duration in milliseconds to human-readable string
 *
 * @param ms - Duration in milliseconds
 * @returns Formatted duration string
 *
 * @example
 * formatDurationMs(1500) // "1.5s"
 * formatDurationMs(90000) // "1m 30s"
 */
export const formatDurationMs = (ms: number): string => {
  return formatDuration(Math.floor(ms / 1000))
}

/**
 * Format currency value
 *
 * @param value - Amount to format
 * @param currency - Currency symbol (default: '$')
 * @param decimals - Number of decimal places (default: 2)
 * @returns Formatted currency string
 *
 * @example
 * formatCurrency(1234.56) // "$1,234.56"
 * formatCurrency(1000, '¥', 0) // "¥1,000"
 */
export const formatCurrency = (
  value: number | null | undefined,
  currency: string = '$',
  decimals: number = 2
): string => {
  if (value === null || value === undefined || isNaN(value)) return '-'
  return `${currency}${formatNumberWithCommas(value, decimals)}`
}

/**
 * Format phone number (basic US format)
 *
 * @param phone - Phone number string
 * @returns Formatted phone number
 *
 * @example
 * formatPhoneNumber('1234567890') // "(123) 456-7890"
 */
export const formatPhoneNumber = (phone: string | null | undefined): string => {
  if (!phone) return '-'

  const cleaned = phone.replace(/\D/g, '')

  if (cleaned.length === 10) {
    return `(${cleaned.slice(0, 3)}) ${cleaned.slice(3, 6)}-${cleaned.slice(6)}`
  }

  if (cleaned.length === 11 && cleaned[0] === '1') {
    return `+1 (${cleaned.slice(1, 4)}) ${cleaned.slice(4, 7)}-${cleaned.slice(7)}`
  }

  return phone
}

/**
 * Truncate string with ellipsis
 *
 * @param str - String to truncate
 * @param maxLength - Maximum length before truncation
 * @param suffix - Suffix to add when truncated (default: '...')
 * @returns Truncated string
 *
 * @example
 * truncate('This is a long string', 10) // "This is a..."
 * truncate('Short', 10) // "Short"
 */
export const truncate = (
  str: string | null | undefined,
  maxLength: number,
  suffix: string = '...'
): string => {
  if (!str) return ''
  if (str.length <= maxLength) return str
  return str.substring(0, maxLength - suffix.length) + suffix
}

/**
 * Format boolean to Yes/No string
 *
 * @param value - Boolean value
 * @param labels - Custom labels for true/false (default: ['Yes', 'No'])
 * @returns Formatted string
 *
 * @example
 * formatBoolean(true) // "Yes"
 * formatBoolean(false, ['Active', 'Inactive']) // "Inactive"
 */
export const formatBoolean = (
  value: boolean | null | undefined,
  labels: [string, string] = ['Yes', 'No']
): string => {
  if (value === null || value === undefined) return '-'
  return value ? labels[0] : labels[1]
}

/**
 * Format array to comma-separated string
 *
 * @param arr - Array to format
 * @param max - Maximum items to show before truncation (default: all)
 * @returns Comma-separated string
 *
 * @example
 * formatArray(['apple', 'banana', 'cherry']) // "apple, banana, cherry"
 * formatArray(['a', 'b', 'c', 'd'], 2) // "a, b, +2 more"
 */
export const formatArray = (arr: any[] | null | undefined, max?: number): string => {
  if (!arr || arr.length === 0) return '-'

  if (max && arr.length > max) {
    const visible = arr.slice(0, max)
    const remaining = arr.length - max
    return `${visible.join(', ')}, +${remaining} more`
  }

  return arr.join(', ')
}

/**
 * Format bytes per second to human-readable transfer rate
 *
 * @param bytesPerSecond - Transfer rate in bytes per second
 * @returns Formatted transfer rate string
 *
 * @example
 * formatTransferRate(1048576) // "1.00 MB/s"
 * formatTransferRate(5120) // "5.00 KB/s"
 */
export const formatTransferRate = (bytesPerSecond: number): string => {
  return `${formatFileSize(bytesPerSecond)}/s`
}

/**
 * Format storage quota display (used / total)
 *
 * @param used - Used storage in bytes
 * @param total - Total storage in bytes
 * @returns Formatted string like "50.0 GB / 100.0 GB (50%)"
 *
 * @example
 * formatStorageQuota(50 * 1024 ** 3, 100 * 1024 ** 3) // "50.00 GB / 100.00 GB (50.0%)"
 */
export const formatStorageQuota = (used: number, total: number): string => {
  const usedStr = formatFileSize(used)
  const totalStr = formatFileSize(total)
  const percent = formatPercent(used / total)
  return `${usedStr} / ${totalStr} (${percent})`
}
