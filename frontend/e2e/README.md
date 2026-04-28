# E2E Tests (Playwright)

End-to-end tests that drive a real Chromium browser against the rsbuild
dev server. We picked Playwright because:

- Multi-browser, headless, parallel — fastest E2E option in 2026
- 1st-class TypeScript, no transpiler shenanigans
- Auto-waiting + auto-retry built in (no `await page.waitForTimeout(500)`)
- **Codegen**: record clicks → generate test code, see "Recording new
  tests" below

## Running

```bash
# Run all tests headless. Auto-starts rsbuild dev on port 4173.
npm run test:e2e

# Open Playwright's time-travel debugger (recommended for writing/fixing)
npm run test:e2e:ui

# Watch the browser do the test (slower, useful for visual debug)
npm run test:e2e:headed

# Open the last HTML report (with traces + screenshots)
npm run test:e2e:report
```

First-time setup (CI does this automatically):

```bash
npm run test:e2e:install   # downloads Chromium (~150MB) into ~/.cache/ms-playwright
```

## Recording new tests (no manual selector hunting)

This is the killer feature. Start codegen, click your way through a flow,
get TypeScript code you can paste into a new `*.spec.ts`:

```bash
# Start the dev server in one terminal:
npx rsbuild dev --port 4173

# In another terminal, start codegen — it opens a new browser window
# alongside the inspector:
npm run test:e2e:codegen
```

Codegen automatically prefers stable selectors (`getByRole`, `getByLabel`,
`getByPlaceholder`) over CSS, so generated code stays readable. Save the
emitted code to `e2e/<flow-name>.spec.ts` and you're done.

## Mocking the backend

Most tests should NOT depend on a running backend — that's slow and
flaky. Use `page.route()` to intercept network requests and return
canned JSON. See `login.spec.ts` for the pattern: stub
`POST /api/v1/auth/login` and `GET /api/v1/users/me` to drive the
client through a successful-login state without postgres or fastapi.

For tests that genuinely need a backend (rare — API-contract tests
should live in the backend repo), set `E2E_BASE_URL` to point at a
running stack:

```bash
E2E_BASE_URL=http://localhost:3000 npm run test:e2e
```

## Why port 4173?

rsbuild's default dev port (3000) and Vite's (5173) are commonly taken
by other apps on developer machines. We use 4173 to avoid the clash.
Override via `E2E_PORT=<n>` if 4173 is also taken.

## Why `NO_PROXY=localhost`?

Some dev machines have `http_proxy` env vars set globally (clash, v2ray,
corporate proxies). Without `NO_PROXY=localhost,127.0.0.1`, axios and
Playwright route localhost calls through the proxy, which often
returns 503 because it can't reach the dev server. The npm scripts in
`package.json` set this for you.
