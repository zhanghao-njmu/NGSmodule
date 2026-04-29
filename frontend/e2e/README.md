# E2E Tests (Playwright + Console-Error Sentinel)

End-to-end tests that drive a real Chromium browser against the rsbuild
dev server, hitting a **real backend** stack (postgres + redis + minio +
fastapi via `docker-compose.smoketest.yml`).

## Why this stack

Per the `frontend-testing-stack` skill (origin: a real bug-discovery
session) the high-leverage rule is:

> Mock the backend in Playwright = anti-pattern. The whole reason
> Playwright catches bugs Vitest can't is because it uses the **real
> backend** — that's how it surfaces stale endpoints, deploy
> mismatches, role-gated 403s, schema drift, auth refresh races.

Each spec attaches a **console-error sentinel** that fails loud on:

- 4xx/5xx responses from `/api/v1/*` (excluding 401, which is
  intentional auth probing)
- Console errors (excluding rsbuild HMR + antd shim warnings)
- Uncaught JS exceptions (`pageerror`)

This caught **4 real production bugs** the day it was introduced:

1. `dist/index.html` had a Vite-leftover `<script src="/src/main.tsx">`
   that would 404 in prod (dev server's SPA fallback masked it).
2. `dashboard /api/v1/stats/quick` returned 500 with
   `type object 'File' has no attribute 'size'` — `File.size` was
   wrong, the SQLAlchemy column is `file_size` (4 callsites in
   stats_service.py + analytics_service.py).
3. Backend `/auth/register` rate limits at 3/min — without a cached
   E2E user the suite tripped the limiter on the second spec.
4. `dist/index.html` hard-coded `<script src="/src/main.tsx">` —
   browser MIME mismatch caught by sentinel on a static-page spec.

## Pre-flight (one time per machine)

```bash
npm run test:e2e:install   # download Chromium (~150MB) → ~/.cache/ms-playwright
```

## Running locally

```bash
# 1. Start the smoke docker stack (publishes backend on 127.0.0.1:8765).
cd ..  # → repo root
docker compose -p ngs-smoke -f docker-compose.yml \
  -f docker-compose.smoketest.yml up -d

# Apply migrations once after a fresh `down -v`:
docker exec ngs-smoke-backend alembic upgrade head

# 2. Run the suite. Playwright auto-starts rsbuild dev on port 4173
#    with BACKEND_PROXY_TARGET=http://127.0.0.1:8765 so the browser
#    runs same-origin and never trips backend CORS.
cd frontend
npm run test:e2e            # headless, all specs
npm run test:e2e:ui         # time-travel debugger
npm run test:e2e:headed     # watch the browser
npm run test:e2e:report     # open last HTML report
```

Tear down when done:

```bash
docker compose -p ngs-smoke -f docker-compose.yml \
  -f docker-compose.smoketest.yml down -v
```

## Recording new tests with Codegen

Click your way through a flow → get TypeScript code:

```bash
# In one terminal, start the dev server pointing at the smoke backend:
BACKEND_PROXY_TARGET=http://127.0.0.1:8765 \
  VITE_API_URL=/api/v1 \
  npx rsbuild dev --port 4173

# In another, start codegen:
npm run test:e2e:codegen
```

Codegen prefers stable selectors (`getByRole`, `getByPlaceholder`).
Paste the output into a new `e2e/<flow-name>.spec.ts`. Wrap each spec
with the sentinel — see `console-errors.spec.ts` for the pattern.

## Anti-patterns the sentinel refuses

- **Mocking the backend** with `page.route()` — turns Playwright
  into slow Vitest and silences the bugs it exists to catch. If a
  test is flaky from backend slowness, fix the slowness, don't mock.
- **Static-navigation specs** — `goto + waitForLoadState + assert
sentinel.all() == []` passes against a stuck Loading screen.
  Trigger real interactions (click, fill form, navigate menus) so
  the fetch chain runs. Verify your spec by intentionally breaking
  an endpoint once — sentinel must turn red.
- **Adding 4xx other than 401 to the ignore list** — 403/404/500
  are signal. Only filter dev-only noise (HMR) and deliberate
  auth probes (401).

## Why `http_proxy=` empty in npm scripts

Many lab/dev boxes have `http_proxy=127.0.0.1:10809` (clash/v2ray)
in their global env. Chromium subprocess inherits it and tries to
"proxy" localhost traffic, which fails (Network Error / 503). Two
defenses combined:

1. `http_proxy= https_proxy= HTTP_PROXY= HTTPS_PROXY=` in npm
   scripts unsets them for Playwright + spawned processes.
2. `--proxy-server=direct://` + `--proxy-bypass-list=*` in
   `playwright.config.ts` `launchOptions.args` hard-disables proxy
   resolution at the Chromium layer.

## Why port 4173

rsbuild's default 3000 and Vite's 5173 are commonly squatted on dev
boxes. 4173 (Vite preview default) was clear. Override via
`E2E_PORT=<n>` if you need to.
