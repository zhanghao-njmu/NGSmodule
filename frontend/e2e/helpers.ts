import { test as base, expect, type Page, type ConsoleMessage } from '@playwright/test'

/**
 * E2E helpers for NGSmodule.
 *
 *  attachErrorSentinel(page) — auto-fail any spec where the backend
 *    returns 4xx/5xx, the page logs console errors, or there's an
 *    uncaught JS exception. The high-leverage piece of the testing
 *    stack: surfaces deploy/role/contract bugs that the spec isn't
 *    directly asserting against.
 *
 *  registerAndLogin(page) — register one E2E user (cached across the
 *    test run because backend rate-limits register at 3/min), log in,
 *    and plant the JWT into the frontend origin's localStorage using
 *    the EXACT keys src/utils/tokenManager.ts reads. Off-by-one keys
 *    = silent auth failure with no console error.
 *
 * Subtleties
 * ----------
 * - Frontend talks to backend through the rsbuild dev proxy (set by
 *   playwright.config.ts → rsbuild reads BACKEND_PROXY_TARGET) so the
 *   browser runs same-origin and never trips backend CORS.
 * - tokenManager validates JWT structure (atob payload). Random-string
 *   tokens fail with "Invalid token format" — we use the real backend's
 *   real JWT.
 * - localStorage is per-origin. We page.goto a same-origin route BEFORE
 *   localStorage.setItem, then page.reload() so App.tsx's checkAuth
 *   useEffect runs again on the fresh mount and sees the planted token.
 */

export interface ErrorSentinel {
  consoleErrors: string[]
  pageErrors: string[]
  badResponses: string[]
  summary(): string
  all(): string[]
}

export function attachErrorSentinel(
  page: Page,
  options: {
    ignoreUrls?: RegExp[]
    ignoreStatuses?: number[]
    ignoreConsolePatterns?: RegExp[]
  } = {},
): ErrorSentinel {
  const consoleErrors: string[] = []
  const pageErrors: string[] = []
  const badResponses: string[] = []

  const ignoreUrls = [
    // rsbuild HMR socket — fails in headless / not part of prod.
    /\/rsbuild-hmr/,
    /\/sockjs-node/,
    /\/__webpack_hmr/,
    ...(options.ignoreUrls || []),
  ]
  // 401 is intentional auth probing in apiClient. DO NOT add 403/404/500
  // here — those are signal.
  const ignoreStatuses = options.ignoreStatuses || [401]
  const ignoreConsolePatterns = [
    /\[Fast Refresh\]/,
    /webpack-hmr/i,
    /HMR/,
    /\[antd:/,
    ...(options.ignoreConsolePatterns || []),
  ]

  page.on('console', (msg: ConsoleMessage) => {
    if (msg.type() !== 'error') return
    const text = msg.text()
    if (ignoreConsolePatterns.some((re) => re.test(text))) return
    consoleErrors.push(text)
  })
  page.on('pageerror', (err) => {
    pageErrors.push(err.message)
  })
  page.on('response', (resp) => {
    const status = resp.status()
    if (status < 400) return
    if (ignoreStatuses.includes(status)) return
    if (ignoreUrls.some((re) => re.test(resp.url()))) return
    badResponses.push(`${status} ${resp.request().method()} ${resp.url()}`)
  })

  return {
    consoleErrors,
    pageErrors,
    badResponses,
    all() {
      return [...consoleErrors, ...pageErrors, ...badResponses]
    },
    summary() {
      const lines: string[] = []
      if (consoleErrors.length) {
        lines.push('CONSOLE ERRORS:')
        lines.push(...consoleErrors.map((e) => '  ' + e))
      }
      if (pageErrors.length) {
        lines.push('PAGE ERRORS:')
        lines.push(...pageErrors.map((e) => '  ' + e))
      }
      if (badResponses.length) {
        lines.push('BAD RESPONSES (4xx/5xx, 401 + HMR ignored):')
        lines.push(...badResponses.map((r) => '  ' + r))
      }
      return lines.join('\n') || '(no errors)'
    },
  }
}

const API_BASE = process.env.E2E_API_BASE_URL || 'http://127.0.0.1:8765'
const REGISTER_PATH = '/api/v1/auth/register'
const LOGIN_PATH = '/api/v1/auth/login'

// Match src/utils/tokenManager.ts.
const TOKEN_KEY = 'auth_token'
const TOKEN_EXPIRY_KEY = 'auth_token_expiry'

// Cache one E2E user across the whole test run because backend
// rate-limits /auth/register at 3/min. Sequential workers + this cache
// = 1 register call per `npx playwright test` run.
let cachedUser: { username: string; email: string; password: string; token: string } | null = null

/**
 * Register + login one shared E2E user (cached). Plant the JWT into
 * the frontend origin's localStorage and reload so AuthGuard remounts
 * with an authed store.
 */
export async function registerAndLogin(
  page: Page,
  options: { forceFresh?: boolean; prefix?: string } = {},
): Promise<{ username: string; email: string; password: string; token: string }> {
  const request = page.context().request

  if (options.forceFresh || !cachedUser) {
    const prefix = options.prefix || 'e2e'
    const stamp = `${Date.now()}${Math.floor(Math.random() * 1000)}`
    const username = `${prefix}${stamp}`
    const email = `${prefix}${stamp}@example.com`
    // Backend password validator wants ≥12 chars + mixed classes.
    const password = 'E2eTest12345!'

    // Step 1: register. 400/409 (duplicate) and 429 (rate limit) are
    // tolerated — we'll fall through to login as long as a user exists.
    const reg = await request.post(`${API_BASE}${REGISTER_PATH}`, {
      data: { username, email, password },
    })
    if (!reg.ok() && reg.status() !== 400 && reg.status() !== 409 && reg.status() !== 429) {
      throw new Error(`E2E register failed: ${reg.status()} ${await reg.text()}`)
    }

    // Step 2: login (form-data, FastAPI OAuth2 password flow).
    const formBody = `username=${encodeURIComponent(username)}&password=${encodeURIComponent(password)}`
    const login = await request.post(`${API_BASE}${LOGIN_PATH}`, {
      data: formBody,
      headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
    })
    if (!login.ok()) {
      throw new Error(`E2E login failed: ${login.status()} ${await login.text()}`)
    }
    const { access_token } = (await login.json()) as { access_token: string }
    cachedUser = { username, email, password, token: access_token }
  }

  const creds = cachedUser

  // Step 3: plant on the frontend origin (localStorage is per-origin).
  await page.goto('/login')
  await page.evaluate(
    ([token, k, ek]) => {
      localStorage.setItem(k, token)
      try {
        const payload = JSON.parse(atob(token.split('.')[1]))
        if (payload.exp) localStorage.setItem(ek, String(payload.exp * 1000))
      } catch {
        /* tokenManager will decode on demand */
      }
    },
    [creds.token, TOKEN_KEY, TOKEN_EXPIRY_KEY],
  )

  // App.tsx's checkAuth() useEffect already fired on the first mount
  // (before token was planted). Reload remounts with the token in place.
  await page.reload()

  return creds
}

export const test = base
export { expect }
