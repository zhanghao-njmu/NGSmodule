/**
 * Rsbuild configuration.
 *
 * Drop-in Vite replacement for the NGSmodule SPA. Rsbuild is built on
 * Rspack (Rust-powered Webpack) and produces ~5-10x faster builds while
 * keeping the dev experience close to Vite (lazy compilation + ESM HMR).
 *
 * Migration notes:
 *   - `loadEnv({ prefixes: ['VITE_', 'PUBLIC_'] })` lets us keep existing
 *     VITE_* env vars without renaming.
 *   - `import.meta.env.DEV` / `.PROD` / `.MODE` are exposed natively.
 *   - `@/...` alias mirrors the previous Vite + tsconfig paths setup.
 *   - `manualChunks` map ports verbatim from the old vite.config.ts.
 */
import { defineConfig, loadEnv } from '@rsbuild/core'
import { pluginReact } from '@rsbuild/plugin-react'
import path from 'path'

const { publicVars } = loadEnv({ prefixes: ['VITE_', 'PUBLIC_'] })

export default defineConfig({
  plugins: [pluginReact()],

  source: {
    entry: { index: './src/main.tsx' },
    // Expose VITE_* / PUBLIC_* env vars on import.meta.env
    define: publicVars,
  },

  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },

  html: {
    template: './index.html',
  },

  server: {
    port: 3000,
    host: '0.0.0.0',
    proxy: {
      // Default to localhost:8000 for `npm run dev`. E2E sets
      // BACKEND_PROXY_TARGET to point at the smoke docker stack
      // (typically http://127.0.0.1:8765) so tests run same-origin
      // and avoid the backend's CORS preflight altogether.
      '/api': {
        target: process.env.BACKEND_PROXY_TARGET || 'http://localhost:8000',
        changeOrigin: true,
      },
      '/api/v1/ws': {
        target: (process.env.BACKEND_PROXY_TARGET || 'http://localhost:8000').replace(/^http/, 'ws'),
        ws: true,
        changeOrigin: true,
      },
    },
  },

  output: {
    distPath: { root: 'dist' },
    sourceMap: { js: false },
    cleanDistPath: true,
  },

  performance: {
    chunkSplit: {
      strategy: 'split-by-experience',
      override: {
        cacheGroups: {
          'react-vendor': {
            test: /[\\/]node_modules[\\/](react|react-dom|react-router-dom|scheduler)[\\/]/,
            name: 'react-vendor',
            chunks: 'all',
            priority: 30,
          },
          'antd-vendor': {
            test: /[\\/]node_modules[\\/](antd|@ant-design)[\\/]/,
            name: 'antd-vendor',
            chunks: 'all',
            priority: 25,
          },
          'chart-vendor': {
            test: /[\\/]node_modules[\\/](plotly\.js|react-plotly\.js)[\\/]/,
            name: 'chart-vendor',
            chunks: 'all',
            priority: 20,
          },
          'query-vendor': {
            test: /[\\/]node_modules[\\/](@tanstack)[\\/]/,
            name: 'query-vendor',
            chunks: 'all',
            priority: 20,
          },
        },
      },
    },
  },
})
