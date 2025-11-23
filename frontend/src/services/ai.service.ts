/**
 * AI Intelligence Service
 * Handles all AI-powered features including recommendations, QC, anomaly detection, and grouping
 */
import apiClient from './api'
import type {
  ParameterRecommendation,
  PipelineRecommendation,
  RecommendationRequest,
  QCReport,
  AutoQCRequest,
  Anomaly,
  AnomalyDetectionRequest,
  AnomalyDetectionReport,
  SmartGroupingRequest,
  SmartGroupingResult,
  AIAssistantMessage,
  AIAssistantConversation,
  AnalysisInsight,
  InsightsReport,
  ResourcePrediction,
  SuccessPrediction,
  AISystemStatus,
} from '@/types/ai'

class AIService {
  // ============================================================
  // Parameter Recommendations
  // ============================================================

  /**
   * Get pipeline parameter recommendations based on context
   */
  async getParameterRecommendations(request: RecommendationRequest): Promise<PipelineRecommendation> {
    const recommendation = await apiClient.post<PipelineRecommendation>('/ai/recommendations/parameters', request)
    return recommendation
  }

  /**
   * Get recommendations for specific parameter
   */
  async getParameterOptions(
    pipelineType: string,
    parameterName: string,
    context?: Record<string, any>,
  ): Promise<ParameterRecommendation> {
    const recommendation = await apiClient.post<ParameterRecommendation>(
      `/ai/recommendations/parameters/${pipelineType}/${parameterName}`,
      { context },
    )
    return recommendation
  }

  /**
   * Get similar successful pipeline runs
   */
  async getSimilarRuns(request: RecommendationRequest): Promise<
    Array<{
      id: string
      similarity: number
      parameters: Record<string, any>
      outcome: string
    }>
  > {
    const runs = await apiClient.post('/ai/recommendations/similar-runs', request)
    return runs
  }

  // ============================================================
  // Quality Control
  // ============================================================

  /**
   * Run automatic quality control analysis
   */
  async runAutoQC(request: AutoQCRequest): Promise<QCReport> {
    const report = await apiClient.post<QCReport>('/ai/qc/auto-analyze', request)
    return report
  }

  /**
   * Get QC recommendations for a sample
   */
  async getQCRecommendations(sampleId: string): Promise<string[]> {
    const recommendations = await apiClient.get<string[]>(`/ai/qc/recommendations/${sampleId}`)
    return recommendations
  }

  /**
   * Batch QC analysis for multiple samples
   */
  async batchQCAnalysis(sampleIds: string[]): Promise<QCReport[]> {
    const reports = await apiClient.post<QCReport[]>('/ai/qc/batch-analyze', {
      sampleIds,
    })
    return reports
  }

  /**
   * Predict QC issues before sequencing
   */
  async predictQCIssues(
    metadata: Record<string, any>,
  ): Promise<Array<{ issue: string; probability: number; prevention: string }>> {
    const predictions = await apiClient.post('/ai/qc/predict-issues', { metadata })
    return predictions
  }

  // ============================================================
  // Anomaly Detection
  // ============================================================

  /**
   * Detect anomalies in samples or project
   */
  async detectAnomalies(request: AnomalyDetectionRequest): Promise<AnomalyDetectionReport> {
    const report = await apiClient.post<AnomalyDetectionReport>('/ai/anomaly/detect', request)
    return report
  }

  /**
   * Get anomaly details
   */
  async getAnomalyDetails(anomalyId: string): Promise<Anomaly> {
    const anomaly = await apiClient.get<Anomaly>(`/ai/anomaly/${anomalyId}`)
    return anomaly
  }

  /**
   * Apply automatic fix for anomaly
   */
  async applyAnomalyFix(anomalyId: string, fixAction: string): Promise<{ success: boolean; message: string }> {
    const result = await apiClient.post(`/ai/anomaly/${anomalyId}/fix`, {
      action: fixAction,
    })
    return result
  }

  /**
   * Monitor for real-time anomalies
   */
  async monitorAnomalies(projectId: string, _callback: (anomaly: Anomaly) => void): Promise<void> {
    // TODO: Implement WebSocket connection for real-time monitoring
    console.log('Monitoring anomalies for project:', projectId)
    // Placeholder for WebSocket implementation
  }

  // ============================================================
  // Smart Sample Grouping
  // ============================================================

  /**
   * Perform smart sample grouping
   */
  async smartGroupSamples(request: SmartGroupingRequest): Promise<SmartGroupingResult> {
    const result = await apiClient.post<SmartGroupingResult>('/ai/grouping/smart-group', request)
    return result
  }

  /**
   * Suggest comparison groups
   */
  async suggestComparisons(projectId: string): Promise<
    Array<{
      group1: string
      group2: string
      comparisonType: string
      reason: string
      confidence: number
    }>
  > {
    const suggestions = await apiClient.get(`/ai/grouping/suggest-comparisons/${projectId}`)
    return suggestions
  }

  /**
   * Validate sample grouping
   */
  async validateGrouping(groups: Array<{ name: string; sampleIds: string[] }>): Promise<{
    valid: boolean
    issues: string[]
    suggestions: string[]
  }> {
    const validation = await apiClient.post('/ai/grouping/validate', { groups })
    return validation
  }

  // ============================================================
  // AI Assistant
  // ============================================================

  /**
   * Send message to AI assistant
   */
  async sendAssistantMessage(
    conversationId: string,
    message: string,
    context?: Record<string, any>,
  ): Promise<AIAssistantMessage> {
    const response = await apiClient.post<AIAssistantMessage>(
      `/ai/assistant/conversations/${conversationId}/messages`,
      { message, context },
    )
    return response
  }

  /**
   * Create new conversation
   */
  async createConversation(title: string, context?: Record<string, any>): Promise<AIAssistantConversation> {
    const conversation = await apiClient.post<AIAssistantConversation>('/ai/assistant/conversations', {
      title,
      context,
    })
    return conversation
  }

  /**
   * Get conversation history
   */
  async getConversation(conversationId: string): Promise<AIAssistantConversation> {
    const conversation = await apiClient.get<AIAssistantConversation>(`/ai/assistant/conversations/${conversationId}`)
    return conversation
  }

  /**
   * List all conversations
   */
  async listConversations(): Promise<AIAssistantConversation[]> {
    const conversations = await apiClient.get<AIAssistantConversation[]>('/ai/assistant/conversations')
    return conversations
  }

  // ============================================================
  // Analysis Insights
  // ============================================================

  /**
   * Get AI-generated insights for a project
   */
  async getProjectInsights(projectId: string): Promise<InsightsReport> {
    const insights = await apiClient.get<InsightsReport>(`/ai/insights/project/${projectId}`)
    return insights
  }

  /**
   * Get insights for specific analysis
   */
  async getAnalysisInsights(analysisType: string, data: any): Promise<AnalysisInsight[]> {
    const insights = await apiClient.post<AnalysisInsight[]>('/ai/insights/analyze', {
      type: analysisType,
      data,
    })
    return insights
  }

  /**
   * Mark insight as reviewed
   */
  async markInsightReviewed(insightId: string): Promise<void> {
    await apiClient.put(`/ai/insights/${insightId}/reviewed`)
  }

  // ============================================================
  // Predictions
  // ============================================================

  /**
   * Predict resource requirements for pipeline
   */
  async predictResources(
    pipelineType: string,
    sampleCount: number,
    dataSize: number,
    parameters?: Record<string, any>,
  ): Promise<ResourcePrediction> {
    const prediction = await apiClient.post<ResourcePrediction>('/ai/predictions/resources', {
      pipelineType,
      sampleCount,
      dataSize,
      parameters,
    })
    return prediction
  }

  /**
   * Predict success probability for analysis
   */
  async predictSuccess(analysisConfig: Record<string, any>): Promise<SuccessPrediction> {
    const prediction = await apiClient.post<SuccessPrediction>('/ai/predictions/success', analysisConfig)
    return prediction
  }

  /**
   * Predict optimal analysis timeline
   */
  async predictTimeline(projectId: string): Promise<{
    totalDuration: string
    milestones: Array<{ name: string; date: string; confidence: number }>
  }> {
    const timeline = await apiClient.get(`/ai/predictions/timeline/${projectId}`)
    return timeline
  }

  // ============================================================
  // System Status
  // ============================================================

  /**
   * Get AI system status and capabilities
   */
  async getSystemStatus(): Promise<AISystemStatus> {
    const status = await apiClient.get<AISystemStatus>('/ai/status')
    return status
  }

  /**
   * Check if AI features are available
   */
  async isAvailable(): Promise<boolean> {
    try {
      const status = await this.getSystemStatus()
      return status.status === 'operational'
    } catch (error) {
      console.error('Failed to check AI availability:', error)
      return false
    }
  }

  // ============================================================
  // Feedback & Learning
  // ============================================================

  /**
   * Submit feedback on AI recommendation
   */
  async submitFeedback(
    recommendationId: string,
    feedback: {
      helpful: boolean
      accuracy?: number
      comments?: string
    },
  ): Promise<void> {
    await apiClient.post(`/ai/feedback/${recommendationId}`, feedback)
  }

  /**
   * Report incorrect prediction
   */
  async reportIncorrectPrediction(predictionId: string, actualOutcome: any, comments?: string): Promise<void> {
    await apiClient.post(`/ai/feedback/prediction/${predictionId}`, {
      actualOutcome,
      comments,
    })
  }
}

export const aiService = new AIService()
export default aiService
