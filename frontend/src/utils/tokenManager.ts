/**
 * Token Manager - Secure token storage and management
 *
 * Security improvements:
 * 1. Centralized token management
 * 2. Token validation before use
 * 3. Automatic cleanup on expiration
 * 4. Event-based token change notifications
 *
 * Note: For production, consider using httpOnly cookies with backend support
 * This implementation uses localStorage with additional security measures
 */

const TOKEN_KEY = 'auth_token'
const TOKEN_EXPIRY_KEY = 'auth_token_expiry'

// Event for token changes
const tokenChangeEvent = new CustomEvent('tokenChange')

/**
 * Parse JWT token payload without verification
 * Only used for checking expiration client-side
 */
function parseJwtPayload(token: string): { exp?: number } | null {
  try {
    const parts = token.split('.')
    if (parts.length !== 3) {
      return null
    }
    const payload = JSON.parse(atob(parts[1]))
    return payload
  } catch {
    return null
  }
}

/**
 * Check if token is expired
 */
function isTokenExpired(token: string): boolean {
  const payload = parseJwtPayload(token)
  if (!payload?.exp) {
    return true
  }

  // Add 60 second buffer for clock skew
  const expiryTime = payload.exp * 1000
  return Date.now() >= expiryTime - 60000
}

/**
 * Token Manager singleton
 */
export const tokenManager = {
  /**
   * Get the stored token
   * Returns null if token is expired or invalid
   */
  getToken(): string | null {
    try {
      const token = localStorage.getItem(TOKEN_KEY)
      if (!token) {
        return null
      }

      // Check if token is expired
      if (isTokenExpired(token)) {
        this.clearToken()
        return null
      }

      return token
    } catch {
      return null
    }
  },

  /**
   * Store a new token
   */
  setToken(token: string): void {
    try {
      // Validate token format
      if (!token || typeof token !== 'string') {
        throw new Error('Invalid token')
      }

      const payload = parseJwtPayload(token)
      if (!payload) {
        throw new Error('Invalid token format')
      }

      localStorage.setItem(TOKEN_KEY, token)

      // Store expiry time for reference
      if (payload.exp) {
        localStorage.setItem(TOKEN_EXPIRY_KEY, String(payload.exp * 1000))
      }

      // Dispatch token change event
      window.dispatchEvent(tokenChangeEvent)
    } catch (error) {
      console.error('Error storing token:', error)
      throw error
    }
  },

  /**
   * Clear stored token
   */
  clearToken(): void {
    try {
      localStorage.removeItem(TOKEN_KEY)
      localStorage.removeItem(TOKEN_EXPIRY_KEY)
      window.dispatchEvent(tokenChangeEvent)
    } catch {
      // Ignore errors during cleanup
    }
  },

  /**
   * Check if user has a valid token
   */
  hasValidToken(): boolean {
    return this.getToken() !== null
  },

  /**
   * Get token expiry time in milliseconds
   */
  getTokenExpiry(): number | null {
    try {
      const expiry = localStorage.getItem(TOKEN_EXPIRY_KEY)
      return expiry ? parseInt(expiry, 10) : null
    } catch {
      return null
    }
  },

  /**
   * Get remaining token validity in milliseconds
   */
  getTokenRemainingTime(): number {
    const expiry = this.getTokenExpiry()
    if (!expiry) {
      return 0
    }
    return Math.max(0, expiry - Date.now())
  },

  /**
   * Subscribe to token changes
   */
  onTokenChange(callback: () => void): () => void {
    const handler = () => callback()
    window.addEventListener('tokenChange', handler)
    return () => window.removeEventListener('tokenChange', handler)
  },
}

export default tokenManager
