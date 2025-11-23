/**
 * AI and Quality Control types
 */

export type QCStatus = 'excellent' | 'good' | 'acceptable' | 'poor' | 'failed'
export type IssueSeverity = 'critical' | 'warning' | 'info'

export interface QCMetric {
  name: string
  value: number
  unit?: string
  status: QCStatus
  explanation: string
  threshold?: {
    excellent: number
    good: number
    acceptable: number
  }
}

export interface QCIssue {
  severity: IssueSeverity
  message: string
  suggestion: string
  affectedMetrics?: string[]
}

export interface QCReport {
  sampleId: string
  sampleName: string
  overallScore: number
  overallStatus: QCStatus
  passedChecks: number
  totalChecks: number
  timestamp: string
  issues: QCIssue[]
  metrics: QCMetric[]
  recommendations?: string[]
  visualizations?: {
    type: string
    data: unknown
    config?: Record<string, unknown>
  }[]
}

export interface AutoQCRequest {
  sampleId: string
  fastqFiles?: string[]
  bamFile?: string
  includeVisualizations?: boolean
}

export interface ParameterRecommendation {
  parameter: string
  recommendedValue: unknown
  confidence: number
  reasoning: string
  reason?: string
  basedOn?: {
    samples: number
    successRate: number
    averageQuality: number
  }
  alternatives?: {
    value: unknown
    confidence: number
    reason: string
  }[]
}

export interface PipelineRecommendation {
  pipelineType: string
  sampleType?: string
  organism?: string
  parameters: ParameterRecommendation[]
  confidence: number
  overallConfidence?: number
  warnings?: string[]
  suggestions?: string[]
  basedOnSamples?: number
  estimatedRuntime?: string
  timestamp: string
}

export interface AIAnalysisResult {
  id: string
  type: 'qc' | 'parameter_optimization' | 'error_diagnosis' | 'performance_tuning'
  status: 'pending' | 'running' | 'completed' | 'failed'
  input: Record<string, unknown>
  output?: Record<string, unknown>
  confidence?: number
  recommendations?: string[]
  timestamp: string
  error?: string
}

export interface AIInsight {
  id: string
  category: 'performance' | 'quality' | 'cost' | 'efficiency'
  title: string
  description: string
  impact: 'high' | 'medium' | 'low'
  actionable: boolean
  action?: {
    label: string
    callback: () => void
  }
  timestamp: string
}
