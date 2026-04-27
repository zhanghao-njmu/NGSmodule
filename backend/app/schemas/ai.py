"""
AI service schemas for intelligent features
"""
from pydantic import BaseModel, Field
from typing import List, Dict, Optional, Any
from datetime import datetime
from enum import Enum
import uuid


# ============================================================================
# Enums
# ============================================================================

class AISystemStatusEnum(str, Enum):
    """AI system operational status"""
    OPERATIONAL = "operational"
    DEGRADED = "degraded"
    MAINTENANCE = "maintenance"
    OFFLINE = "offline"


class QCStatusEnum(str, Enum):
    """QC analysis status"""
    PASS = "pass"
    WARNING = "warning"
    FAIL = "fail"


class AnomalySeverity(str, Enum):
    """Anomaly severity level"""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


class AnomalyStatus(str, Enum):
    """Anomaly status"""
    DETECTED = "detected"
    INVESTIGATING = "investigating"
    RESOLVED = "resolved"
    DISMISSED = "dismissed"


class InsightCategory(str, Enum):
    """Insight category"""
    QUALITY = "quality"
    PERFORMANCE = "performance"
    OPTIMIZATION = "optimization"
    ANOMALY = "anomaly"
    PREDICTION = "prediction"


# ============================================================================
# Parameter Recommendations
# ============================================================================

class RecommendationRequest(BaseModel):
    """Pipeline recommendation request"""
    pipeline_type: str
    sample_count: Optional[int] = None
    data_size: Optional[int] = None
    organism: Optional[str] = None
    sequencing_type: Optional[str] = None
    quality_threshold: Optional[float] = None
    context: Optional[Dict[str, Any]] = None


class ParameterValue(BaseModel):
    """Single parameter recommendation value"""
    name: str
    recommended_value: Any
    confidence: float = Field(..., ge=0.0, le=1.0)
    reasoning: str
    alternative_values: Optional[List[Any]] = None


class ParameterRecommendation(BaseModel):
    """Parameter recommendation for specific parameter"""
    parameter_name: str
    pipeline_type: str
    recommended_value: Any
    confidence: float
    reasoning: str
    alternatives: List[Dict[str, Any]] = []
    based_on_runs: int = 0


class PipelineRecommendation(BaseModel):
    """Complete pipeline parameter recommendations"""
    pipeline_type: str
    parameters: List[ParameterValue]
    overall_confidence: float
    similar_runs_analyzed: int
    estimated_success_rate: float
    estimated_runtime: Optional[float] = None  # in seconds
    notes: List[str] = []
    generated_at: datetime = Field(default_factory=datetime.utcnow)


class SimilarRun(BaseModel):
    """Similar successful pipeline run"""
    id: str
    similarity: float
    parameters: Dict[str, Any]
    outcome: str
    success_metric: Optional[float] = None


# ============================================================================
# Quality Control
# ============================================================================

class AutoQCRequest(BaseModel):
    """Auto QC analysis request"""
    sample_id: str
    file_id: Optional[str] = None
    qc_type: str = "comprehensive"  # basic, comprehensive, custom
    parameters: Optional[Dict[str, Any]] = None


class QCMetric(BaseModel):
    """Single QC metric"""
    name: str
    value: float
    threshold: Optional[float] = None
    status: QCStatusEnum
    description: Optional[str] = None


class QCRecommendation(BaseModel):
    """QC improvement recommendation"""
    issue: str
    severity: str
    recommendation: str
    expected_impact: str


class QCReport(BaseModel):
    """Quality control report"""
    sample_id: str
    overall_status: QCStatusEnum
    overall_score: float = Field(..., ge=0.0, le=100.0)
    metrics: List[QCMetric]
    recommendations: List[QCRecommendation] = []
    summary: str
    generated_at: datetime = Field(default_factory=datetime.utcnow)


class BatchQCRequest(BaseModel):
    """Batch QC request"""
    sampleIds: List[str] = Field(..., alias="sampleIds", min_items=1)

    class Config:
        allow_population_by_field_name = True
        populate_by_name = True


class QCIssuePrediction(BaseModel):
    """Predicted QC issue"""
    issue: str
    probability: float = Field(..., ge=0.0, le=1.0)
    prevention: str
    impact: str


class PredictQCRequest(BaseModel):
    """Predict QC issues request"""
    metadata: Dict[str, Any]


# ============================================================================
# Anomaly Detection
# ============================================================================

class AnomalyDetectionRequest(BaseModel):
    """Anomaly detection request"""
    project_id: Optional[str] = None
    sample_ids: Optional[List[str]] = None
    detection_type: str = "comprehensive"
    sensitivity: str = "medium"  # low, medium, high
    parameters: Optional[Dict[str, Any]] = None


class Anomaly(BaseModel):
    """Detected anomaly"""
    id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    type: str
    severity: AnomalySeverity
    status: AnomalyStatus = AnomalyStatus.DETECTED
    title: str
    description: str
    affected_resource: Dict[str, str]  # {"type": "sample", "id": "..."}
    detected_at: datetime = Field(default_factory=datetime.utcnow)
    metric_values: Optional[Dict[str, Any]] = None
    suggested_actions: List[str] = []
    auto_fixable: bool = False
    confidence: float = Field(..., ge=0.0, le=1.0)


class AnomalyDetectionReport(BaseModel):
    """Anomaly detection report"""
    detection_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    project_id: Optional[str] = None
    anomalies: List[Anomaly]
    total_detected: int
    by_severity: Dict[str, int]
    summary: str
    detected_at: datetime = Field(default_factory=datetime.utcnow)


class AnomalyFixRequest(BaseModel):
    """Anomaly fix request"""
    action: str
    parameters: Optional[Dict[str, Any]] = None


class AnomalyFixResponse(BaseModel):
    """Anomaly fix result"""
    success: bool
    message: str
    anomaly_id: str
    fix_applied: Optional[str] = None
    timestamp: datetime = Field(default_factory=datetime.utcnow)


# ============================================================================
# Smart Sample Grouping
# ============================================================================

class SmartGroupingRequest(BaseModel):
    """Smart grouping request"""
    project_id: str
    sample_ids: Optional[List[str]] = None
    grouping_criteria: List[str] = []  # e.g., ["condition", "tissue", "time_point"]
    min_group_size: int = 3
    max_groups: int = 10


class SampleGroup(BaseModel):
    """Sample group"""
    name: str
    sample_ids: List[str]
    common_attributes: Dict[str, Any]
    confidence: float


class SmartGroupingResult(BaseModel):
    """Smart grouping result"""
    groups: List[SampleGroup]
    ungrouped_samples: List[str] = []
    grouping_quality: float
    total_samples: int
    suggestions: List[str] = []
    generated_at: datetime = Field(default_factory=datetime.utcnow)


class ComparisonSuggestion(BaseModel):
    """Comparison suggestion"""
    group1: str
    group2: str
    comparisonType: str
    reason: str
    confidence: float


class GroupValidationRequest(BaseModel):
    """Group validation request"""
    groups: List[Dict[str, Any]]  # [{"name": "...", "sampleIds": [...]}]


class GroupValidationResult(BaseModel):
    """Group validation result"""
    valid: bool
    issues: List[str] = []
    suggestions: List[str] = []


# ============================================================================
# AI Assistant
# ============================================================================

class AssistantMessageRequest(BaseModel):
    """Assistant message request"""
    message: str
    context: Optional[Dict[str, Any]] = None


class AIAssistantMessage(BaseModel):
    """AI assistant message"""
    id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    conversation_id: str
    role: str  # user, assistant, system
    content: str
    timestamp: datetime = Field(default_factory=datetime.utcnow)
    metadata: Optional[Dict[str, Any]] = None


class AIAssistantConversation(BaseModel):
    """AI assistant conversation"""
    id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    title: str
    messages: List[AIAssistantMessage] = []
    context: Optional[Dict[str, Any]] = None
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)


class CreateConversationRequest(BaseModel):
    """Create conversation request"""
    title: str
    context: Optional[Dict[str, Any]] = None


# ============================================================================
# Analysis Insights
# ============================================================================

class AnalysisInsight(BaseModel):
    """Single analysis insight"""
    id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    category: InsightCategory
    title: str
    description: str
    importance: str  # low, medium, high
    confidence: float
    actions: List[str] = []
    related_data: Optional[Dict[str, Any]] = None
    reviewed: bool = False
    generated_at: datetime = Field(default_factory=datetime.utcnow)


class InsightsReport(BaseModel):
    """Insights report for a project"""
    project_id: str
    insights: List[AnalysisInsight]
    summary: str
    total_insights: int
    by_category: Dict[str, int]
    generated_at: datetime = Field(default_factory=datetime.utcnow)


class AnalyzeInsightsRequest(BaseModel):
    """Insights analysis request"""
    type: str
    data: Dict[str, Any]


# ============================================================================
# Predictions
# ============================================================================

class ResourcePredictionRequest(BaseModel):
    """Resource prediction request"""
    pipelineType: str
    sampleCount: int
    dataSize: int
    parameters: Optional[Dict[str, Any]] = None

    class Config:
        populate_by_name = True


class ResourcePrediction(BaseModel):
    """Resource requirement prediction"""
    cpu_cores: int
    memory_gb: float
    storage_gb: float
    estimated_duration_seconds: float
    confidence: float
    notes: List[str] = []


class SuccessPrediction(BaseModel):
    """Success probability prediction"""
    success_probability: float = Field(..., ge=0.0, le=1.0)
    risk_factors: List[str] = []
    success_factors: List[str] = []
    recommendations: List[str] = []
    confidence: float


class TimelineMilestone(BaseModel):
    """Timeline milestone"""
    name: str
    date: str
    confidence: float


class TimelinePrediction(BaseModel):
    """Project timeline prediction"""
    totalDuration: str
    milestones: List[TimelineMilestone]
    confidence: float
    generated_at: datetime = Field(default_factory=datetime.utcnow)


# ============================================================================
# System Status
# ============================================================================

class AICapability(BaseModel):
    """AI capability information"""
    name: str
    available: bool
    version: Optional[str] = None
    description: Optional[str] = None


class AISystemStatus(BaseModel):
    """AI system status"""
    status: AISystemStatusEnum
    version: str
    capabilities: List[AICapability]
    models_loaded: List[str] = []
    last_updated: datetime = Field(default_factory=datetime.utcnow)
    uptime_seconds: float = 0.0
    notes: List[str] = []


# ============================================================================
# Feedback
# ============================================================================

class FeedbackRequest(BaseModel):
    """Feedback request"""
    helpful: bool
    accuracy: Optional[float] = Field(None, ge=0.0, le=1.0)
    comments: Optional[str] = None


class IncorrectPredictionRequest(BaseModel):
    """Incorrect prediction report"""
    actualOutcome: Any
    comments: Optional[str] = None


class FeedbackResponse(BaseModel):
    """Feedback acknowledgement"""
    success: bool
    message: str
    feedback_id: str = Field(default_factory=lambda: str(uuid.uuid4()))
    timestamp: datetime = Field(default_factory=datetime.utcnow)
