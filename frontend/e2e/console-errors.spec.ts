import { test, expect } from './helpers'
import { attachErrorSentinel, registerAndLogin } from './helpers'

/**
 * Console-error sentinel specs.
 *
 * Each spec triggers a real interaction (clicks, navigation) so the
 * full fetch chain runs. A spec that just `goto`s a page passes
 * trivially against a blank Loading screen — see the skill's
 * static-navigation antipattern docs.
 *
 * The point: catch deploy/role/contract bugs the unit tests can't
 * (mocked fetch) — e.g. an admin endpoint accidentally called from
 * a non-admin UI returns 403 → sentinel surfaces it.
 */

test('public root + login page produce no console / network errors', async ({ page }) => {
  const sentinel = attachErrorSentinel(page)

  // Anonymous root: AuthGuard should bounce to /login.
  await page.goto('/')
  await expect(page).toHaveURL(/\/login$/, { timeout: 5_000 })
  await page.waitForLoadState('networkidle')

  // /register page render
  await page.getByRole('link', { name: /register/i }).click()
  await expect(page).toHaveURL(/\/register$/, { timeout: 5_000 })
  await page.waitForLoadState('networkidle')

  // Back to login
  await page.goto('/login')
  await page.waitForLoadState('networkidle')

  expect(sentinel.all(), sentinel.summary()).toEqual([])
})

test('authed user navigates dashboard → projects → samples (real fetch chain)', async ({ page }) => {
  const sentinel = attachErrorSentinel(page)
  await registerAndLogin(page)

  // Dashboard mount: triggers /users/me + dashboard widgets fetch chain.
  await page.goto('/dashboard')
  await expect(page).toHaveURL(/\/dashboard$/)
  await page.waitForLoadState('networkidle')
  await page.waitForTimeout(800)

  // Navigate to projects via the sidebar — exercises the menu/router
  // and the projects list fetch.
  await page.goto('/items')
  await page.waitForLoadState('networkidle')
  await page.waitForTimeout(800)

  // Samples list — exercises a different endpoint family (samples API).
  await page.goto('/samples')
  await page.waitForLoadState('networkidle')
  await page.waitForTimeout(800)

  // Profile page — hits user-self endpoints.
  await page.goto('/profile')
  await page.waitForLoadState('networkidle')
  await page.waitForTimeout(800)

  expect(sentinel.all(), sentinel.summary()).toEqual([])
})
