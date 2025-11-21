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
        <Route path="/dashboard" element={<Dashboard />} />
        <Route path="/projects" element={<ProjectList />} />
        <Route path="/samples" element={<SampleList />} />
        <Route path="/files" element={<FileList />} />
        <Route path="/tasks" element={<TaskList />} />
      </Route>
    </Routes>
  )
}

export default App
