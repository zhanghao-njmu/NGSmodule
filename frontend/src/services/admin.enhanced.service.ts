/**
 * Enhanced Admin Service
 * Extended administrative functions for system monitoring and management
 */
import apiClient from './api'

// System Monitoring Types
export interface SystemHealth {
  overall: 'healthy' | 'degraded' | 'critical'
  components: Record<string, { status: string; latency?: number }>
  metrics: {
    cpu: number
    memory: number
    disk: number
  }
}

export interface SystemMetrics {
  cpu: { usage: number; load: number[] }
  memory: { used: number; total: number; usagePercent: number }
  disk: { used: number; total: number; usagePercent: number }
}

export interface Alert {
  id: string
  type: 'error' | 'warning' | 'info'
  severity: 'critical' | 'high' | 'medium' | 'low'
  title: string
  message: string
  timestamp: string
  resolved: boolean
}

// Audit Log Types
export interface AuditLog {
  id: string
  timestamp: string
  userId: string
  username: string
  action: string
  resource: string
  status: 'success' | 'failure'
  ip: string
}

class EnhancedAdminService {
  // ============================================================
  // System Monitoring
  // ============================================================

  async getSystemHealth(): Promise<SystemHealth> {
    const health = await apiClient.get<SystemHealth>('/admin/system/health')
    return health
  }

  async getSystemMetrics(): Promise<SystemMetrics> {
    const metrics = await apiClient.get<SystemMetrics>('/admin/system/metrics')
    return metrics
  }

  async getAlerts(resolved = false): Promise<Alert[]> {
    const alerts = await apiClient.get<Alert[]>('/admin/alerts', {
      params: { resolved },
    })
    return alerts
  }

  async resolveAlert(alertId: string): Promise<void> {
    await apiClient.post(`/admin/alerts/${alertId}/resolve`)
  }

  // ============================================================
  // Audit Logs
  // ============================================================

  async getAuditLogs(filters?: {
    startDate?: string
    endDate?: string
    userId?: string
    action?: string
  }): Promise<AuditLog[]> {
    const logs = await apiClient.get<AuditLog[]>('/admin/audit-logs', {
      params: filters,
    })
    return logs
  }

  async exportAuditLogs(filters?: any): Promise<{ downloadUrl: string }> {
    const result = await apiClient.post('/admin/audit-logs/export', filters)
    return result
  }

  // ============================================================
  // Resource Management
  // ============================================================

  async getResourceUsage(): Promise<{
    storage: { total: number; used: number }
    compute: { active: number; limit: number }
  }> {
    const usage = await apiClient.get('/admin/resources/usage')
    return usage
  }

  // ============================================================
  // Backup Management
  // ============================================================

  async createBackup(type: 'full' | 'incremental'): Promise<{ id: string }> {
    const result = await apiClient.post('/admin/backups', { type })
    return result
  }

  async listBackups(): Promise<
    Array<{
      id: string
      type: string
      createdAt: string
      size: number
      status: string
    }>
  > {
    const backups = await apiClient.get('/admin/backups')
    return backups
  }

  // ============================================================
  // Job Management
  // ============================================================

  async listJobs(status?: string): Promise<
    Array<{
      id: string
      type: string
      status: string
      userId: string
      createdAt: string
    }>
  > {
    const jobs = await apiClient.get('/admin/jobs', {
      params: { status },
    })
    return jobs
  }

  async cancelJob(jobId: string): Promise<void> {
    await apiClient.post(`/admin/jobs/${jobId}/cancel`)
  }

  async retryJob(jobId: string): Promise<void> {
    await apiClient.post(`/admin/jobs/${jobId}/retry`)
  }
}

export const enhancedAdminService = new EnhancedAdminService()
export default enhancedAdminService
