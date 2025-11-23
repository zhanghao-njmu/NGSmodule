/**
 * Analytics & Visualization Types
 * Types for advanced analytics, charts, and data visualization
 */

// ============================================================
// Chart Data Types
// ============================================================

export interface ChartData {
  labels: string[]
  datasets: ChartDataset[]
}

export interface ChartDataset {
  label: string
  data: number[]
  backgroundColor?: string | string[]
  borderColor?: string | string[]
  borderWidth?: number
  fill?: boolean
  tension?: number
  pointRadius?: number
  pointHoverRadius?: number
}

export interface ChartOptions {
  responsive?: boolean
  maintainAspectRatio?: boolean
  plugins?: {
    legend?: {
      display?: boolean
      position?: 'top' | 'bottom' | 'left' | 'right'
    }
    title?: {
      display?: boolean
      text?: string
    }
    tooltip?: {
      enabled?: boolean
    }
  }
  scales?: {
    x?: any
    y?: any
  }
}

// ============================================================
// Analytics Dashboard Types
// ============================================================

export interface AnalyticsSummary {
  period: 'today' | 'week' | 'month' | 'year' | 'all'
  projectStats: {
    total: number
    active: number
    completed: number
    failed: number
    successRate: number
  }
  sampleStats: {
    total: number
    processed: number
    pending: number
    avgQualityScore: number
  }
  pipelineStats: {
    totalRuns: number
    successfulRuns: number
    failedRuns: number
    avgRuntime: string
    avgCost: number
  }
  storageStats: {
    totalUsed: number
    totalAvailable: number
    usageByType: Record<string, number>
    trend: 'increasing' | 'decreasing' | 'stable'
  }
  activityStats: {
    dailyActiveUsers: number
    totalSessions: number
    avgSessionDuration: string
    peakUsageHour: number
  }
}

export interface TimeSeriesData {
  timestamp: string
  value: number
  category?: string
  metadata?: Record<string, any>
}

export interface TrendAnalysis {
  metric: string
  currentValue: number
  previousValue: number
  change: number
  changePercent: number
  trend: 'up' | 'down' | 'stable'
  prediction?: {
    nextValue: number
    confidence: number
  }
}

// ============================================================
// Project Comparison Types
// ============================================================

export interface ProjectComparison {
  items: string[]
  metrics: ComparisonMetric[]
  summary: {
    bestPerforming: string
    recommendations: string[]
  }
  visualizations: {
    type: 'bar' | 'line' | 'radar' | 'scatter'
    data: any
    title: string
  }[]
}

export interface ComparisonMetric {
  name: string
  unit?: string
  values: Record<string, number>
  winner?: string
  description: string
}

export interface ProjectPerformance {
  projectId: string
  projectName: string
  scores: {
    qualityScore: number
    efficiencyScore: number
    costEffectivenessScore: number
    timelinessScore: number
    overallScore: number
  }
  rankings: {
    quality: number
    efficiency: number
    cost: number
    timeline: number
    overall: number
  }
  strengths: string[]
  weaknesses: string[]
}

// ============================================================
// Visualization Types
// ============================================================

export type VisualizationType =
  | 'line'
  | 'bar'
  | 'pie'
  | 'doughnut'
  | 'scatter'
  | 'bubble'
  | 'radar'
  | 'polar'
  | 'heatmap'
  | 'boxplot'
  | 'violin'
  | 'sankey'
  | 'treemap'
  | 'network'
  | 'histogram'

export interface Visualization {
  id: string
  type: VisualizationType
  title: string
  description?: string
  data: any
  options?: ChartOptions
  interactive?: boolean
  downloadable?: boolean
  shareable?: boolean
}

export interface HeatmapData {
  rows: string[]
  columns: string[]
  values: number[][]
  colorScale?: {
    min: string
    mid: string
    max: string
  }
}

export interface NetworkData {
  nodes: Array<{
    id: string
    label: string
    size?: number
    color?: string
    group?: string
  }>
  edges: Array<{
    source: string
    target: string
    weight?: number
    label?: string
  }>
}

// ============================================================
// Report Types
// ============================================================

export interface AnalyticsReport {
  id: string
  title: string
  description: string
  type: 'project' | 'sample' | 'pipeline' | 'system' | 'custom'
  generatedAt: string
  generatedBy: string
  period: {
    start: string
    end: string
  }
  sections: ReportSection[]
  attachments?: Array<{
    name: string
    type: string
    url: string
  }>
  exportFormats: ('pdf' | 'excel' | 'csv' | 'json')[]
}

export interface ReportSection {
  title: string
  type: 'text' | 'table' | 'chart' | 'metrics' | 'list'
  content: any
  visualization?: Visualization
  insights?: string[]
}

// ============================================================
// Knowledge Base Types
// ============================================================

export interface KnowledgeArticle {
  id: string
  title: string
  category: KnowledgeCategory
  tags: string[]
  content: string
  summary: string
  author: string
  createdAt: string
  updatedAt: string
  views: number
  helpful: number
  notHelpful: number
  difficulty: 'beginner' | 'intermediate' | 'advanced'
  estimatedReadTime: number
  relatedArticles: string[]
  attachments?: Array<{
    name: string
    type: string
    url: string
  }>
}

export type KnowledgeCategory =
  | 'getting-started'
  | 'pipelines'
  | 'data-management'
  | 'analysis'
  | 'quality-control'
  | 'troubleshooting'
  | 'best-practices'
  | 'api-reference'
  | 'faq'

export interface Tutorial {
  id: string
  title: string
  description: string
  category: KnowledgeCategory
  difficulty: 'beginner' | 'intermediate' | 'advanced'
  estimatedDuration: number
  steps: TutorialStep[]
  prerequisites?: string[]
  outcomes: string[]
  completed: boolean
  progress: number
}

export interface TutorialStep {
  id: string
  title: string
  content: string
  type: 'text' | 'video' | 'interactive' | 'code' | 'quiz'
  media?: {
    type: 'image' | 'video' | 'gif'
    url: string
    caption?: string
  }
  code?: {
    language: string
    content: string
    runnable?: boolean
  }
  quiz?: {
    question: string
    options: string[]
    correctAnswer: number
    explanation: string
  }
  completed: boolean
}

// ============================================================
// Statistics Types
// ============================================================

export interface StatisticalAnalysis {
  metric: string
  sampleSize: number
  mean: number
  median: number
  mode: number
  stdDev: number
  variance: number
  min: number
  max: number
  quartiles: {
    q1: number
    q2: number
    q3: number
  }
  outliers: number[]
  distribution: 'normal' | 'skewed' | 'bimodal' | 'uniform' | 'unknown'
  confidenceInterval: {
    level: number
    lower: number
    upper: number
  }
}

export interface CorrelationMatrix {
  variables: string[]
  correlations: number[][]
  pValues: number[][]
  significanceLevel: number
}

// ============================================================
// Export Types
// ============================================================

export interface ExportRequest {
  type: 'chart' | 'report' | 'data' | 'visualization'
  format: 'png' | 'svg' | 'pdf' | 'excel' | 'csv' | 'json'
  itemId: string
  options?: {
    includeData?: boolean
    includeMetadata?: boolean
    resolution?: 'low' | 'medium' | 'high'
    orientation?: 'portrait' | 'landscape'
  }
}

export interface ExportResult {
  success: boolean
  downloadUrl?: string
  expiresAt?: string
  fileSize?: number
  error?: string
}

// ============================================================
// Dashboard Widget Types
// ============================================================

export interface DashboardWidget {
  id: string
  type: 'chart' | 'metric' | 'table' | 'text' | 'custom'
  title: string
  position: {
    x: number
    y: number
    width: number
    height: number
  }
  config: any
  refreshInterval?: number
  lastUpdated?: string
}

export interface DashboardLayout {
  id: string
  name: string
  userId: string
  widgets: DashboardWidget[]
  isDefault: boolean
  createdAt: string
  updatedAt: string
}
