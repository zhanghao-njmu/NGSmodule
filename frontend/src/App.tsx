import { useEffect } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { ConfigProvider, theme as antdTheme } from 'antd'
import { authStore } from '@/store/authStore'
import { useTheme } from '@/store/themeStore'
import { lightTheme, darkTheme } from '@/styles/theme.config'
import { MainLayout } from '@/layouts/MainLayout'
import { AuthLayout } from '@/layouts/AuthLayout'
import { Login } from '@/pages/auth/Login'
import { Register } from '@/pages/auth/Register'
import { Dashboard } from '@/pages/dashboard/Dashboard'
import { ProjectList } from '@/pages/projects/ProjectList'
import { SampleList } from '@/pages/samples/SampleList'
import { FileList } from '@/pages/files/FileList'
import { TaskList } from '@/pages/tasks/TaskList'
import { PipelineList } from '@/pages/pipelines/PipelineList'
import { ResultDetail } from '@/pages/results/ResultDetail'
import { ResultList } from '@/pages/results/ResultList'
import { AdminDashboard } from '@/pages/admin/AdminDashboard'
import { ProgressBar, ErrorBoundary } from '@/components/common'

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
  const { mode, isDark } = useTheme()

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
            path="/projects"
            element={
              <ErrorBoundary>
                <ProjectList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/samples"
            element={
              <ErrorBoundary>
                <SampleList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/files"
            element={
              <ErrorBoundary>
                <FileList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/pipelines"
            element={
              <ErrorBoundary>
                <PipelineList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/tasks"
            element={
              <ErrorBoundary>
                <TaskList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/results"
            element={
              <ErrorBoundary>
                <ResultList />
              </ErrorBoundary>
            }
          />
          <Route
            path="/results/:id"
            element={
              <ErrorBoundary>
                <ResultDetail />
              </ErrorBoundary>
            }
          />
          <Route
            path="/admin"
            element={
              <ErrorBoundary>
                <AdminDashboard />
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
