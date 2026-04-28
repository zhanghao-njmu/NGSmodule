import { useEffect, lazy, Suspense } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { ConfigProvider, theme as antdTheme, Spin } from 'antd'
import { authStore } from '@/store/authStore'
import { useTheme } from '@/store/themeStore'
import { lightTheme, darkTheme } from '@/styles/theme.config'
import { MainLayout } from '@/layouts/MainLayout'
import { AuthLayout } from '@/layouts/AuthLayout'
import { ProgressBar, ErrorBoundary } from '@/components/common'

// Eager load: Auth pages (needed for initial load)
import { Login } from '@/pages/auth/Login'
import { Register } from '@/pages/auth/Register'
import { Dashboard } from '@/pages/dashboard/Dashboard'

// Lazy load: All other pages (code splitting for better performance)
const ProjectList = lazy(() => import('@/pages/projects/ProjectList').then((m) => ({ default: m.ProjectList })))
const SampleList = lazy(() => import('@/pages/samples/SampleList').then((m) => ({ default: m.SampleList })))
const FileList = lazy(() => import('@/pages/files/FileList').then((m) => ({ default: m.FileList })))
const TaskList = lazy(() => import('@/pages/tasks/TaskList').then((m) => ({ default: m.TaskList })))
const PipelineList = lazy(() => import('@/pages/pipelines/PipelineList').then((m) => ({ default: m.PipelineList })))
const ResultDetail = lazy(() => import('@/pages/results/ResultDetail').then((m) => ({ default: m.ResultDetail })))
const ResultList = lazy(() => import('@/pages/results/ResultList').then((m) => ({ default: m.ResultList })))
const AdminDashboard = lazy(() => import('@/pages/admin/AdminDashboard').then((m) => ({ default: m.AdminDashboard })))
const ProfilePage = lazy(() => import('@/pages/profile/ProfilePage').then((m) => ({ default: m.ProfilePage })))
const SettingsPage = lazy(() => import('@/pages/settings/SettingsPage').then((m) => ({ default: m.SettingsPage })))
const NotificationsPage = lazy(() =>
  import('@/pages/notifications/NotificationsPage').then((m) => ({ default: m.NotificationsPage })),
)
const AIDashboard = lazy(() => import('@/pages/ai/AIDashboard').then((m) => ({ default: m.AIDashboard })))
const AnalyticsDashboard = lazy(() =>
  import('@/pages/analytics/AnalyticsDashboard').then((m) => ({ default: m.AnalyticsDashboard })),
)
const KnowledgeBase = lazy(() => import('@/pages/knowledge/KnowledgeBase').then((m) => ({ default: m.KnowledgeBase })))
const DataDownloadsPage = lazy(() =>
  import('@/pages/data-downloads/DataDownloadsPage').then((m) => ({ default: m.DataDownloadsPage })),
)

// Loading fallback component
const PageLoader = () => (
  <div style={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: '400px' }}>
    <Spin size="large" tip="Loading..." />
  </div>
)

// Protected Route Component
const ProtectedRoute = ({ children }: { children: React.ReactNode }) => {
  const { isAuthenticated } = authStore()

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />
  }

  return <>{children}</>
}

function App() {
  const { checkAuth } = authStore()
  const { isDark } = useTheme()

  useEffect(() => {
    checkAuth()
  }, [checkAuth])

  // Select theme based on mode
  const currentTheme = isDark ? darkTheme : lightTheme

  // Add dark mode algorithm if in dark mode
  const themeConfig = {
    ...currentTheme,
    algorithm: isDark ? antdTheme.darkAlgorithm : antdTheme.defaultAlgorithm,
  }

  return (
    <ConfigProvider theme={themeConfig}>
      <ErrorBoundary onReset={() => window.location.reload()}>
        <ProgressBar />
        <Routes>
          {/* Public Routes */}
          <Route element={<AuthLayout />}>
            <Route path="/login" element={<Login />} />
            <Route path="/register" element={<Register />} />
          </Route>

          {/* Protected Routes */}
          <Route
            element={
              <ProtectedRoute>
                <MainLayout />
              </ProtectedRoute>
            }
          >
            <Route path="/" element={<Navigate to="/dashboard" replace />} />
            <Route
              path="/dashboard"
              element={
                <ErrorBoundary>
                  <Dashboard />
                </ErrorBoundary>
              }
            />
            <Route
              path="/items"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <ProjectList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/samples"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <SampleList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/files"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <FileList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/pipelines"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <PipelineList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/tasks"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <TaskList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/results"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <ResultList />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/results/:id"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <ResultDetail />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/admin"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <AdminDashboard />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/profile"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <ProfilePage />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/settings"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <SettingsPage />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/notifications"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <NotificationsPage />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/ai"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <AIDashboard />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/analytics"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <AnalyticsDashboard />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/knowledge"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <KnowledgeBase />
                  </Suspense>
                </ErrorBoundary>
              }
            />
            <Route
              path="/data-downloads"
              element={
                <ErrorBoundary>
                  <Suspense fallback={<PageLoader />}>
                    <DataDownloadsPage />
                  </Suspense>
                </ErrorBoundary>
              }
            />
          </Route>
        </Routes>
      </ErrorBoundary>
    </ConfigProvider>
  )
}

export default App
