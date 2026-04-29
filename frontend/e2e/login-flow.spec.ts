import { test, expect } from './helpers'
import { attachErrorSentinel, registerAndLogin } from './helpers'

/**
 * Login flow — drives the real login form against the real backend.
 *
 * Spec order: form rendering → validation → end-to-end login redirect.
 */

test('login page renders form fields', async ({ page }) => {
  const sentinel = attachErrorSentinel(page)
  await page.goto('/login')
  await expect(page.getByPlaceholder('Username or Email')).toBeVisible()
  await expect(page.getByPlaceholder('Password')).toBeVisible()
  await expect(page.getByRole('button', { name: /log in|sign in|登录/i })).toBeVisible()
  expect(sentinel.all(), sentinel.summary()).toEqual([])
})

test('empty submit triggers antd Form validation', async ({ page }) => {
  const sentinel = attachErrorSentinel(page)
  await page.goto('/login')
  await page.getByRole('button', { name: /log in|sign in|登录/i }).click()
  await expect(page.getByText('Please input your username or email!')).toBeVisible()
  await expect(page.getByText('Please input your password!')).toBeVisible()
  expect(sentinel.all(), sentinel.summary()).toEqual([])
})

test('register + login → real backend → land on dashboard', async ({ page }) => {
  const sentinel = attachErrorSentinel(page)
  // Register a fresh user via API + plant the JWT (no UI typing); then
  // navigate to the protected root and verify AuthGuard accepts us.
  const { username } = await registerAndLogin(page)

  // Go to dashboard. AuthGuard reads the token, hits /users/me, and
  // either lets us in or punts to /login.
  await page.goto('/dashboard')
  await expect(page).toHaveURL(/\/dashboard$/, { timeout: 10_000 })

  // Wait for post-mount fetches (tanstack-query prefetches) to settle.
  await page.waitForLoadState('networkidle')
  await page.waitForTimeout(500)

  // The dashboard renders a "Welcome back, <username>!" heading once
  // /users/me resolves — confirms we ARE authenticated, not just
  // URL-matching a public page.
  await expect(page.getByRole('heading', { name: new RegExp(`Welcome back,\\s*${username}`) })).toBeVisible({
    timeout: 10_000,
  })

  expect(sentinel.all(), sentinel.summary()).toEqual([])
})
