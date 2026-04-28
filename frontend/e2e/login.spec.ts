import { test, expect } from '@playwright/test'

/**
 * Smoke E2E: login page + form validation + happy-path login.
 *
 * Notes
 * -----
 * 1. We mock the backend at the network layer (`page.route`) so this test
 *    can run on CI without a postgres/redis/backend stack — just rsbuild
 *    dev server + chromium.
 * 2. The auth service expects:
 *      POST /api/v1/auth/login  -> { access_token, token_type }
 *      GET  /api/v1/auth/me     -> User
 *    Both are intercepted in the happy-path test.
 */

test.describe('login page', () => {
  test('renders form fields and submit button', async ({ page }) => {
    await page.goto('/login')
    await expect(page.getByPlaceholder('Username or Email')).toBeVisible()
    await expect(page.getByPlaceholder('Password')).toBeVisible()
    await expect(page.getByRole('button', { name: /log in|sign in|登录/i })).toBeVisible()
  })

  test('empty submit triggers antd Form validation', async ({ page }) => {
    await page.goto('/login')
    await page.getByRole('button', { name: /log in|sign in|登录/i }).click()
    // antd surfaces the rule message inline
    await expect(page.getByText('Please input your username or email!')).toBeVisible()
    await expect(page.getByText('Please input your password!')).toBeVisible()
  })

  test('successful login redirects to dashboard', async ({ page }) => {
    // The token manager parses access_token client-side (atob the JWT
    // payload), so we need a real JWT shape — a literal random string is
    // rejected as "Invalid token format". Header + payload below decode to
    //   header  = {"alg":"HS256","typ":"JWT"}
    //   payload = {"sub":"00000000-0000-0000-0000-000000000001","exp":9999999999}
    // Signature is irrelevant on the client (the backend would verify it).
    const HEADER = 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9'
    const PAYLOAD = 'eyJzdWIiOiIwMDAwMDAwMC0wMDAwLTAwMDAtMDAwMC0wMDAwMDAwMDAwMDEiLCJleHAiOjk5OTk5OTk5OTl9'
    const FAKE_JWT = `${HEADER}.${PAYLOAD}.signature-not-verified-client-side`

    // Stub backend: return token + a minimal User shape that satisfies the
    // store's expectations (id, username, email, role, is_active).
    await page.route('**/api/v1/auth/login', async (route) => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          access_token: FAKE_JWT,
          token_type: 'bearer',
        }),
      })
    })

    // The auth service calls GET /users/me (not /auth/me) after login.
    await page.route('**/api/v1/users/me', async (route) => {
      await route.fulfill({
        status: 200,
        contentType: 'application/json',
        body: JSON.stringify({
          id: '00000000-0000-0000-0000-000000000001',
          username: 'e2e-user',
          email: 'e2e@test.local',
          full_name: 'E2E Tester',
          role: 'user',
          is_active: true,
          storage_quota: 107374182400,
          storage_used: 0,
          created_at: '2026-01-01T00:00:00',
        }),
      })
    })

    await page.goto('/login')
    await page.getByPlaceholder('Username or Email').fill('e2e-user')
    await page.getByPlaceholder('Password').fill('whatever-stubbed')
    await page.getByRole('button', { name: /log in|sign in|登录/i }).click()

    await expect(page).toHaveURL(/\/dashboard$/, { timeout: 10_000 })
  })
})
