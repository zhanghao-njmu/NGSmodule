/**
 * Admin types and interfaces
 */

export type SystemStatus = 'healthy' | 'degraded' | 'critical' | 'operational'
export type AlertSeverity = 'critical' | 'high' | 'medium' | 'low'
export type AlertType = 'error' | 'warning' | 'info'

export interface ComponentHealth {
  status: SystemStatus
  latency?: number
  message?: string
}

export interface SystemHealth {
  overall: SystemStatus
  components: {
    api: ComponentHealth
    database: ComponentHealth
    storage: ComponentHealth
    queue: ComponentHealth
    [key: string]: ComponentHealth
  }
  metrics: {
    cpu: number
    memory: number
    disk: number
  }
}

export interface SystemMetrics {
  cpu: {
    usage: number
    load: number[]
  }
  memory: {
    used: number
    total: number
    usagePercent: number
  }
  disk: {
    used: number
    total: number
    usagePercent: number
  }
  network?: {
    bytesIn: number
    bytesOut: number
  }
}

export interface SystemAlert {
  id: string
  type: AlertType
  severity: AlertSeverity
  title: string
  message: string
  timestamp: string
  resolved: boolean
  resolvedAt?: string
  resolvedBy?: string
}

export interface SystemLog {
  id: string
  level: 'debug' | 'info' | 'warn' | 'error'
  message: string
  timestamp: string
  source?: string
  metadata?: Record<string, unknown>
}

export interface BackupStatus {
  id: string
  type: 'full' | 'incremental' | 'differential'
  status: 'pending' | 'running' | 'completed' | 'failed'
  startedAt?: string
  completedAt?: string
  size?: number
  location?: string
  error?: string
}

export interface UserActivity {
  userId: string
  username: string
  action: string
  resource: string
  timestamp: string
  ipAddress?: string
  userAgent?: string
}

export interface User {
  id: string
  username: string
  email: string
  full_name?: string
  role: 'admin' | 'user'
  is_active: boolean
  organization?: string
  storage_used?: number
  storage_quota?: number
  created_at: string
  last_login?: string
}

export interface UserAdminUpdate {
  username?: string
  email?: string
  full_name?: string
  role?: 'admin' | 'user'
  is_active?: boolean
  password?: string
  storage_quota?: number
}

export interface SystemStats {
  total_users: number
  active_users: number
  total_projects: number
  total_samples: number
  total_tasks: number
  running_tasks: number
  completed_tasks?: number
  failed_tasks?: number
  storage_used: number
  storage_total: number
  total_storage_used?: number
  total_storage_quota?: number
  cpu_usage?: number
  memory_usage?: number
}

export interface AdminStats {
  totalUsers: number
  activeUsers: number
  totalProjects: number
  totalSamples: number
  totalTasks: number
  runningTasks: number
  storageUsed: number
  storageTotal: number
}

export interface UserStats {
  userId: string
  username: string
  totalProjects: number
  totalSamples: number
  totalTasks: number
  storageUsed: number
  storageQuota: number
  lastActive?: string
  createdAt: string
}
