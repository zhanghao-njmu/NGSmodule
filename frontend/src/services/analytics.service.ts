/**
 * Analytics Service
 * Handles analytics, visualization, and reporting functionality
 */
import apiClient from './api'
import type {
  AnalyticsSummary,
  TimeSeriesData,
  TrendAnalysis,
  ProjectComparison,
  ProjectPerformance,
  Visualization,
  AnalyticsReport,
  KnowledgeArticle,
  KnowledgeCategory,
  Tutorial,
  StatisticalAnalysis,
  CorrelationMatrix,
  ExportRequest,
  ExportResult,
  DashboardLayout,
} from '@/types/analytics'

class AnalyticsService {
  // ============================================================
  // Dashboard Analytics
  // ============================================================

  /**
   * Get analytics summary for dashboard
   */
  async getAnalyticsSummary(
    period: 'today' | 'week' | 'month' | 'year' | 'all' = 'week'
  ): Promise<AnalyticsSummary> {
    const summary = await apiClient.get<AnalyticsSummary>('/analytics/summary', {
      params: { period },
    })
    return summary
  }

  /**
   * Get time series data for metrics
   */
  async getTimeSeriesData(
    metric: string,
    startDate: string,
    endDate: string,
    interval: 'hour' | 'day' | 'week' | 'month' = 'day'
  ): Promise<TimeSeriesData[]> {
    const data = await apiClient.get<TimeSeriesData[]>('/analytics/timeseries', {
      params: { metric, startDate, endDate, interval },
    })
    return data
  }

  /**
   * Get trend analysis for metrics
   */
  async getTrendAnalysis(metrics: string[]): Promise<TrendAnalysis[]> {
    const trends = await apiClient.post<TrendAnalysis[]>('/analytics/trends', {
      metrics,
    })
    return trends
  }

  // ============================================================
  // Project Analytics
  // ============================================================

  /**
   * Compare multiple projects
   */
  async compareProjects(projectIds: string[]): Promise<ProjectComparison> {
    const comparison = await apiClient.post<ProjectComparison>(
      '/analytics/projects/compare',
      { projectIds }
    )
    return comparison
  }

  /**
   * Get project performance analysis
   */
  async getProjectPerformance(projectId: string): Promise<ProjectPerformance> {
    const performance = await apiClient.get<ProjectPerformance>(
      `/analytics/projects/${projectId}/performance`
    )
    return performance
  }

  /**
   * Get project statistics
   */
  async getProjectStatistics(
    projectId: string
  ): Promise<{
    samples: StatisticalAnalysis
    quality: StatisticalAnalysis
    runtime: StatisticalAnalysis
  }> {
    const stats = await apiClient.get(`/analytics/projects/${projectId}/statistics`)
    return stats
  }

  // ============================================================
  // Visualizations
  // ============================================================

  /**
   * Get visualization by ID
   */
  async getVisualization(visualizationId: string): Promise<Visualization> {
    const viz = await apiClient.get<Visualization>(
      `/analytics/visualizations/${visualizationId}`
    )
    return viz
  }

  /**
   * Create custom visualization
   */
  async createVisualization(config: {
    type: string
    title: string
    dataSource: string
    options: any
  }): Promise<Visualization> {
    const viz = await apiClient.post<Visualization>(
      '/analytics/visualizations',
      config
    )
    return viz
  }

  /**
   * Get correlation matrix
   */
  async getCorrelationMatrix(
    projectId: string,
    variables: string[]
  ): Promise<CorrelationMatrix> {
    const matrix = await apiClient.post<CorrelationMatrix>(
      `/analytics/projects/${projectId}/correlation`,
      { variables }
    )
    return matrix
  }

  // ============================================================
  // Reports
  // ============================================================

  /**
   * Generate analytics report
   */
  async generateReport(config: {
    type: string
    projectId?: string
    period: { start: string; end: string }
    sections: string[]
  }): Promise<AnalyticsReport> {
    const report = await apiClient.post<AnalyticsReport>('/analytics/reports', config)
    return report
  }

  /**
   * Get report by ID
   */
  async getReport(reportId: string): Promise<AnalyticsReport> {
    const report = await apiClient.get<AnalyticsReport>(
      `/analytics/reports/${reportId}`
    )
    return report
  }

  /**
   * List all reports
   */
  async listReports(
    filters?: {
      type?: string
      startDate?: string
      endDate?: string
    }
  ): Promise<AnalyticsReport[]> {
    const reports = await apiClient.get<AnalyticsReport[]>('/analytics/reports', {
      params: filters,
    })
    return reports
  }

  /**
   * Export report
   */
  async exportReport(
    reportId: string,
    format: 'pdf' | 'excel' | 'csv'
  ): Promise<ExportResult> {
    const result = await apiClient.post<ExportResult>(
      `/analytics/reports/${reportId}/export`,
      { format }
    )
    return result
  }

  // ============================================================
  // Knowledge Base
  // ============================================================

  /**
   * Search knowledge base
   */
  async searchKnowledge(
    query: string,
    category?: KnowledgeCategory
  ): Promise<KnowledgeArticle[]> {
    const articles = await apiClient.get<KnowledgeArticle[]>('/knowledge/search', {
      params: { q: query, category },
    })
    return articles
  }

  /**
   * Get article by ID
   */
  async getArticle(articleId: string): Promise<KnowledgeArticle> {
    const article = await apiClient.get<KnowledgeArticle>(
      `/knowledge/articles/${articleId}`
    )
    return article
  }

  /**
   * Get articles by category
   */
  async getArticlesByCategory(
    category: KnowledgeCategory
  ): Promise<KnowledgeArticle[]> {
    const articles = await apiClient.get<KnowledgeArticle[]>(
      `/knowledge/categories/${category}`
    )
    return articles
  }

  /**
   * Get popular articles
   */
  async getPopularArticles(limit = 10): Promise<KnowledgeArticle[]> {
    const articles = await apiClient.get<KnowledgeArticle[]>(
      '/knowledge/popular',
      {
        params: { limit },
      }
    )
    return articles
  }

  /**
   * Mark article as helpful
   */
  async markArticleHelpful(articleId: string, helpful: boolean): Promise<void> {
    await apiClient.post(`/knowledge/articles/${articleId}/feedback`, {
      helpful,
    })
  }

  // ============================================================
  // Tutorials
  // ============================================================

  /**
   * Get all tutorials
   */
  async getTutorials(category?: KnowledgeCategory): Promise<Tutorial[]> {
    const tutorials = await apiClient.get<Tutorial[]>('/knowledge/tutorials', {
      params: { category },
    })
    return tutorials
  }

  /**
   * Get tutorial by ID
   */
  async getTutorial(tutorialId: string): Promise<Tutorial> {
    const tutorial = await apiClient.get<Tutorial>(
      `/knowledge/tutorials/${tutorialId}`
    )
    return tutorial
  }

  /**
   * Mark tutorial step as completed
   */
  async completeTutorialStep(
    tutorialId: string,
    stepId: string
  ): Promise<void> {
    await apiClient.post(
      `/knowledge/tutorials/${tutorialId}/steps/${stepId}/complete`
    )
  }

  /**
   * Get user's tutorial progress
   */
  async getTutorialProgress(): Promise<
    Array<{ tutorialId: string; progress: number; completed: boolean }>
  > {
    const progress = await apiClient.get('/knowledge/tutorials/progress')
    return progress
  }

  // ============================================================
  // Data Export
  // ============================================================

  /**
   * Export chart or visualization
   */
  async exportVisualization(request: ExportRequest): Promise<ExportResult> {
    const result = await apiClient.post<ExportResult>(
      '/analytics/export',
      request
    )
    return result
  }

  /**
   * Export data to CSV
   */
  async exportDataToCSV(
    dataType: string,
    itemId: string,
    options?: any
  ): Promise<ExportResult> {
    const result = await apiClient.post<ExportResult>('/analytics/export/csv', {
      dataType,
      itemId,
      options,
    })
    return result
  }

  /**
   * Export multiple items in batch
   */
  async batchExport(requests: ExportRequest[]): Promise<ExportResult[]> {
    const results = await apiClient.post<ExportResult[]>(
      '/analytics/export/batch',
      { requests }
    )
    return results
  }

  // ============================================================
  // Dashboard Customization
  // ============================================================

  /**
   * Get user's dashboard layouts
   */
  async getDashboardLayouts(): Promise<DashboardLayout[]> {
    const layouts = await apiClient.get<DashboardLayout[]>('/analytics/dashboards')
    return layouts
  }

  /**
   * Save dashboard layout
   */
  async saveDashboardLayout(layout: Partial<DashboardLayout>): Promise<DashboardLayout> {
    const saved = await apiClient.post<DashboardLayout>(
      '/analytics/dashboards',
      layout
    )
    return saved
  }

  /**
   * Update dashboard layout
   */
  async updateDashboardLayout(
    layoutId: string,
    updates: Partial<DashboardLayout>
  ): Promise<DashboardLayout> {
    const updated = await apiClient.put<DashboardLayout>(
      `/analytics/dashboards/${layoutId}`,
      updates
    )
    return updated
  }

  /**
   * Delete dashboard layout
   */
  async deleteDashboardLayout(layoutId: string): Promise<void> {
    await apiClient.delete(`/analytics/dashboards/${layoutId}`)
  }

  /**
   * Set default dashboard
   */
  async setDefaultDashboard(layoutId: string): Promise<void> {
    await apiClient.post(`/analytics/dashboards/${layoutId}/default`)
  }
}

export const analyticsService = new AnalyticsService()
export default analyticsService
