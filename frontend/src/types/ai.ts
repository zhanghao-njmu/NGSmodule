/**
 * AI Intelligence Types
 * Types for AI-powered features including recommendations, QC, and analysis
 */

// ============================================================
// Parameter Recommendation Types
// ============================================================

export interface ParameterRecommendation {
  parameter: string
  recommendedValue: any
  confidence: number // 0-1
  reason: string
  alternatives?: Array<{
    value: any
    confidence: number
    reason: string
  }>
  basedOn: {
    similarProjects?: number
    historicalData?: number
    publicDatasets?: number
    bestPractices?: boolean
  }
}

export interface PipelineRecommendation {
  pipelineType: string
  pipelineName: string
  parameters: ParameterRecommendation[]
  overallConfidence: number
  estimatedRuntime: string
  estimatedCost?: number
  warnings?: string[]
  suggestions?: string[]
}

export interface RecommendationRequest {
  pipelineType: string
  sampleType?: string
  sequencingPlatform?: string
  organism?: string
  analysisGoal?: string
  dataQuality?: {
    avgQuality: number
    totalReads: number
    readLength: number
  }
  customContext?: Record<string, any>
}

// ============================================================
// Quality Control Types
// ============================================================

export type QCStatus = 'excellent' | 'good' | 'acceptable' | 'poor' | 'failed'

export interface QCMetric {
  name: string
  value: number
  unit?: string
  threshold: {
    excellent: number
    good: number
    acceptable: number
  }
  status: QCStatus
  explanation: string
}

export interface QCReport {
  sampleId: string
  sampleName: string
  overallScore: number // 0-100
  overallStatus: QCStatus
  timestamp: string
  metrics: QCMetric[]
  recommendations: string[]
  issues: Array<{
    severity: 'critical' | 'warning' | 'info'
    message: string
    suggestion?: string
  }>
  passedChecks: number
  totalChecks: number
}

export interface AutoQCRequest {
  sampleId: string
  fastqFiles?: string[]
  bamFile?: string
  includeVisualizations?: boolean
}

// ============================================================
// Anomaly Detection Types
// ============================================================

export type AnomalyType =
  | 'quality_drop'
  | 'contamination'
  | 'batch_effect'
  | 'outlier_expression'
  | 'sequencing_bias'
  | 'coverage_issue'
  | 'adapter_content'
  | 'gc_bias'

export type AnomalySeverity = 'critical' | 'high' | 'medium' | 'low'

export interface Anomaly {
  id: string
  type: AnomalyType
  severity: AnomalySeverity
  confidence: number // 0-1
  affectedSamples: string[]
  description: string
  detectedAt: string
  metrics: Record<string, number>
  visualization?: {
    type: 'scatter' | 'heatmap' | 'line' | 'bar'
    data: any
  }
  recommendations: string[]
  autoFixAvailable: boolean
  fixActions?: Array<{
    action: string
    description: string
    impact: string
  }>
}

export interface AnomalyDetectionRequest {
  projectId?: string
  sampleIds?: string[]
  checkTypes?: AnomalyType[]
  sensitivity?: 'low' | 'medium' | 'high'
}

export interface AnomalyDetectionReport {
  timestamp: string
  anomalies: Anomaly[]
  summary: {
    total: number
    bySeverity: Record<AnomalySeverity, number>
    byType: Record<AnomalyType, number>
  }
  healthScore: number // 0-100
  recommendations: string[]
}

// ============================================================
// Smart Sample Grouping Types
// ============================================================

export interface SampleGroup {
  id: string
  name: string
  samples: string[]
  groupingReason: string
  confidence: number
  characteristics: Record<string, any>
  suggestedComparisons?: Array<{
    groupId: string
    comparisonType: string
    reason: string
  }>
}

export interface SmartGroupingRequest {
  projectId: string
  samples: Array<{
    id: string
    name: string
    metadata: Record<string, any>
    qcMetrics?: Record<string, number>
  }>
  groupingStrategy?:
    | 'metadata'
    | 'qc_similarity'
    | 'expression_pattern'
    | 'auto'
  maxGroupsPerType?: number
}

export interface SmartGroupingResult {
  groups: SampleGroup[]
  ungroupedSamples: string[]
  groupingStrategy: string
  confidence: number
  visualizations?: Array<{
    type: 'pca' | 'heatmap' | 'dendrogram'
    data: any
    title: string
  }>
  suggestions: string[]
}

// ============================================================
// AI Assistant Types
// ============================================================

export interface AIAssistantMessage {
  id: string
  role: 'user' | 'assistant'
  content: string
  timestamp: string
  context?: {
    projectId?: string
    sampleId?: string
    pipelineId?: string
  }
  suggestions?: string[]
  actions?: Array<{
    label: string
    action: string
    parameters?: Record<string, any>
  }>
}

export interface AIAssistantConversation {
  id: string
  title: string
  messages: AIAssistantMessage[]
  createdAt: string
  updatedAt: string
  context: Record<string, any>
}

// ============================================================
// Analysis Insights Types
// ============================================================

export interface AnalysisInsight {
  id: string
  type:
    | 'differential_expression'
    | 'pathway_enrichment'
    | 'variant_calling'
    | 'quality_issue'
    | 'optimization'
    | 'discovery'
  title: string
  description: string
  confidence: number
  importance: 'critical' | 'high' | 'medium' | 'low'
  affectedSamples?: string[]
  affectedGenes?: string[]
  visualization?: {
    type: string
    data: any
  }
  relatedInsights?: string[]
  actionable: boolean
  suggestedActions?: string[]
  references?: Array<{
    title: string
    url: string
    source: string
  }>
}

export interface InsightsReport {
  projectId: string
  timestamp: string
  insights: AnalysisInsight[]
  summary: {
    totalInsights: number
    byType: Record<string, number>
    byImportance: Record<string, number>
    actionableCount: number
  }
  topRecommendations: string[]
}

// ============================================================
// Prediction Types
// ============================================================

export interface ResourcePrediction {
  estimatedRuntime: {
    value: number
    unit: 'minutes' | 'hours' | 'days'
    confidence: number
  }
  estimatedCost: {
    value: number
    currency: string
    breakdown: Record<string, number>
    confidence: number
  }
  estimatedStorage: {
    value: number
    unit: 'GB' | 'TB'
    confidence: number
  }
  basedOn: {
    similarRuns: number
    dataSize: number
    pipelineComplexity: string
  }
  recommendations: string[]
}

export interface SuccessPrediction {
  overallProbability: number // 0-1
  factors: Array<{
    factor: string
    impact: 'positive' | 'negative' | 'neutral'
    weight: number
    description: string
  }>
  riskFactors: string[]
  recommendations: string[]
}

// ============================================================
// AI Model Info Types
// ============================================================

export interface AIModel {
  name: string
  version: string
  type: 'recommendation' | 'qc' | 'anomaly' | 'grouping' | 'prediction'
  description: string
  trainedOn: {
    datasets: number
    samples: number
    lastUpdated: string
  }
  performance: {
    accuracy?: number
    precision?: number
    recall?: number
    f1Score?: number
  }
  status: 'active' | 'training' | 'deprecated'
}

export interface AISystemStatus {
  available: boolean
  models: AIModel[]
  capabilities: string[]
  usage: {
    requestsToday: number
    averageResponseTime: number
    errorRate: number
  }
  limits: {
    maxRequestsPerDay: number
    maxRequestsPerHour: number
  }
}
