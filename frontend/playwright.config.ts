import { defineConfig, devices } from '@playwright/test'

/**
 * Playwright config for NGSmodule frontend E2E tests.
 *
 * Modes
 * -----
 * Default (`npm run test:e2e`)     auto-starts rsbuild dev. Frontend proxies
 *                                  /api -> http://127.0.0.1:8765 (the smoke
 *                                  docker stack's published backend port).
 * Codegen (`npm run test:e2e:codegen`)
 *                                  opens browser inspector; record clicks
 *                                  -> Playwright code.
 *
 * Required preconditions
 * ----------------------
 * 1. The smoke docker stack is running:
 *      docker compose -p ngs-smoke -f docker-compose.yml \
 *        -f docker-compose.smoketest.yml up -d
 *    (verify: `curl http://127.0.0.1:8765/health` -> {"status":"healthy"})
 * 2. Migrations applied:
 *      docker exec ngs-smoke-backend alembic upgrade head
 *
 * Why these knobs
 * ---------------
 * - `--no-proxy-server` Chromium flag: bypass any inherited http_proxy
 *   env var (clash/v2ray/corp proxy) which silently rewrites localhost
 *   traffic and stalls the dev server. Symptom is a hung Loading screen
 *   with no console errors.
 * - port 4173: rsbuild's default 3000 and Vite's 5173 are commonly
 *   squatted by other apps on dev boxes; 4173 (Vite preview default)
 *   has been clear in our environment.
 * - fullyParallel:false: specs share state (one E2E user registered up
 *   front, multiple specs use it). Sequential keeps the auth user
 *   reusable across specs.
 */
const E2E_PORT = Number(process.env.E2E_PORT || 4173)
const E2E_BASE_URL = process.env.E2E_BASE_URL || `http://127.0.0.1:${E2E_PORT}`

export default defineConfig({
  testDir: './e2e',
  timeout: 30_000,
  expect: { timeout: 5_000 },
  fullyParallel: false,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 1 : 0,
  workers: 1,
  reporter: process.env.CI ? 'list' : [['list'], ['html', { open: 'never' }]],
  use: {
    baseURL: E2E_BASE_URL,
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
    video: 'retain-on-failure',
  },
  projects: [
    {
      name: 'chromium',
      use: {
        ...devices['Desktop Chrome'],
        // The shell's http_proxy=127.0.0.1:10809 (clash/v2ray) gets
        // inherited by Chromium even with NO_PROXY set; it tries to
        // proxy localhost traffic, axios sees ECONNREFUSED-as-503 and
        // the page sits on a Loading screen forever. The chromium-
        // specific flag --proxy-server=direct:// hard-disables proxy
        // resolution at the browser layer.
        launchOptions: {
          args: ['--proxy-server=direct://', '--proxy-bypass-list=*'],
        },
      },
    },
  ],
  webServer: process.env.E2E_BASE_URL
    ? undefined
    : {
        // Run frontend same-origin — rsbuild dev proxies /api/v1/* to
        // the backend, so the browser never makes a cross-origin call
        // and we don't have to mess with BACKEND_CORS_ORIGINS for tests.
        // - VITE_API_URL set to bare /api/v1 → axios uses relative path
        //   → request hits 4173 (frontend dev server) → proxied to 8765.
        // - BACKEND_PROXY_TARGET is read by rsbuild.config.ts to pick
        //   the proxy upstream.
        command: `npx rsbuild dev --port ${E2E_PORT}`,
        url: E2E_BASE_URL,
        reuseExistingServer: false,
        timeout: 120_000,
        env: {
          BACKEND_PROXY_TARGET: 'http://127.0.0.1:8765',
          VITE_API_URL: '/api/v1',
        },
      },
})
