import { defineConfig, devices } from '@playwright/test'

/**
 * Playwright config for NGSmodule frontend E2E tests.
 *
 * Local: `npm run test:e2e` — auto-starts the rsbuild dev server.
 * Record: `npm run test:e2e:codegen` — opens a browser; every click you
 *   make is transcribed into TypeScript test code, ready to paste into
 *   a new file in `e2e/`.
 * UI debug: `npm run test:e2e:ui` — Playwright's time-travel debugger.
 *
 * Port choice: 4173 instead of rsbuild's default 3000 because the dev
 * machine often has other Next.js/Vite apps on 3000/5173. The test
 * server is started on its own port via `rsbuild dev --port`, fully
 * isolated from anything the developer might already be running.
 */
const E2E_PORT = Number(process.env.E2E_PORT || 4173)
const E2E_BASE_URL = process.env.E2E_BASE_URL || `http://localhost:${E2E_PORT}`

export default defineConfig({
  testDir: './e2e',
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 2 : 0,
  workers: process.env.CI ? 1 : undefined,
  reporter: [['html', { open: 'never' }], ['list']],
  use: {
    baseURL: E2E_BASE_URL,
    trace: 'on-first-retry',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
  webServer: process.env.E2E_BASE_URL
    ? undefined
    : {
        command: `npx rsbuild dev --port ${E2E_PORT}`,
        url: E2E_BASE_URL,
        reuseExistingServer: false,
        timeout: 120_000,
      },
})
