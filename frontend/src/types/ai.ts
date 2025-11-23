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
    similarProjects?: number
    bestPractices?: boolean
    publicDatasets?: number
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

export interface RecommendationRequest {
  pipelineType: string
  sampleType?: string
  organism?: string
  existingParameters?: Record<string, unknown>
}

export interface Anomaly {
  id: string
  type: 'outlier' | 'pattern' | 'drift' | 'missing'
  severity: 'critical' | 'warning' | 'info'
  description: string
  affectedSamples?: string[]
  suggestedAction?: string
  detectedAt: string
}

export interface AnomalyDetectionRequest {
  projectId?: string
  sampleIds?: string[]
  metrics?: string[]
  sensitivity?: 'low' | 'medium' | 'high'
}

export interface AnomalyDetectionReport {
  requestId: string
  projectId?: string
  totalSamples: number
  anomaliesDetected: number
  anomalies: Anomaly[]
  summary: string
  timestamp: string
}

export interface SmartGroupingRequest {
  sampleIds: string[]
  strategy?: 'similarity' | 'metadata' | 'quality' | 'automatic'
  minGroupSize?: number
  maxGroups?: number
}

export interface SmartGroupingResult {
  groups: {
    id: string
    name: string
    samples: string[]
    characteristics: Record<string, unknown>
    confidence: number
  }[]
  strategy: string
  timestamp: string
}

export interface AIAssistantMessage {
  id: string
  role: 'user' | 'assistant' | 'system'
  content: string
  timestamp: string
  metadata?: Record<string, unknown>
}

export interface AIAssistantConversation {
  id: string
  userId: string
  projectId?: string
  messages: AIAssistantMessage[]
  title?: string
  createdAt: string
  updatedAt: string
}

export interface AnalysisInsight {
  id: string
  type: 'optimization' | 'warning' | 'tip' | 'best_practice'
  category: string
  title: string
  description: string
  impact: 'high' | 'medium' | 'low'
  actionable: boolean
  recommendation?: string
  relatedEntities?: string[]
}

export interface InsightsReport {
  projectId?: string
  timestamp: string
  totalInsights: number
  insights: AnalysisInsight[]
  summary: string
  categories: Record<string, number>
}

export interface ResourcePrediction {
  estimatedDuration: number
  estimatedMemory: number
  estimatedCPU: number
  estimatedStorage: number
  confidence: number
  basedOn: {
    samples: number
    similarPipelines: number
  }
}

export interface SuccessPrediction {
  successProbability: number
  failureRisks: {
    factor: string
    probability: number
    mitigation?: string
  }[]
  recommendations: string[]
  confidence: number
}

export interface AISystemStatus {
  status: 'operational' | 'degraded' | 'offline'
  availableFeatures: string[]
  queueSize: number
  averageResponseTime: number
  lastUpdate: string
}
