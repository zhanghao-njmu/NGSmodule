import { useEffect } from 'react'
import { Routes, Route, Navigate } from 'react-router-dom'
import { authStore } from '@/store/authStore'
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

  useEffect(() => {
    checkAuth()
  }, [checkAuth])

  return (
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
  )
}

export default App
