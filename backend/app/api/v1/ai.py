"""
AI Intelligence API Router
Handles AI-powered features including recommendations, QC, anomaly detection, and grouping
"""
from typing import List, Dict, Any
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from app.core.deps import get_db, get_current_user
from app.models.user import User
from app.services.ai_service import ai_service
from app.schemas.ai import (
    # Parameter Recommendations
    RecommendationRequest, ParameterRecommendation, PipelineRecommendation, SimilarRun,
    # Quality Control
    AutoQCRequest, QCReport, BatchQCRequest, QCIssuePrediction, PredictQCRequest,
    # Anomaly Detection
    AnomalyDetectionRequest, Anomaly, AnomalyDetectionReport,
    AnomalyFixRequest, AnomalyFixResponse,
    # Smart Sample Grouping
    SmartGroupingRequest, SmartGroupingResult, ComparisonSuggestion,
    GroupValidationRequest, GroupValidationResult,
    # AI Assistant
    AssistantMessageRequest, AIAssistantMessage, AIAssistantConversation,
    CreateConversationRequest,
    # Analysis Insights
    AnalysisInsight, InsightsReport, AnalyzeInsightsRequest,
    # Predictions
    ResourcePredictionRequest, ResourcePrediction, SuccessPrediction, TimelinePrediction,
    # System Status
    AISystemStatus,
    # Feedback
    FeedbackRequest, IncorrectPredictionRequest, FeedbackResponse,
)

router = APIRouter()


# ============================================================================
# Parameter Recommendations
# ============================================================================

@router.post("/recommendations/parameters", response_model=PipelineRecommendation)
async def get_parameter_recommendations(
    request: RecommendationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get pipeline parameter recommendations based on context

    Analyzes similar successful runs and provides AI-powered parameter suggestions.
    """
    return ai_service.get_parameter_recommendations(request, db)


@router.post(
    "/recommendations/parameters/{pipeline_type}/{parameter_name}",
    response_model=ParameterRecommendation
)
async def get_parameter_options(
    pipeline_type: str,
    parameter_name: str,
    request: Dict[str, Any] = None,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get recommendations for specific parameter

    Provides detailed options and alternatives for a specific pipeline parameter.
    """
    context = request.get("context") if request else None
    return ai_service.get_parameter_options(pipeline_type, parameter_name, context, db)


@router.post("/recommendations/similar-runs", response_model=List[SimilarRun])
async def get_similar_runs(
    request: RecommendationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get similar successful pipeline runs

    Finds and returns pipeline runs with similar characteristics and successful outcomes.
    """
    return ai_service.get_similar_runs(request, db)


# ============================================================================
# Quality Control
# ============================================================================

@router.post("/qc/auto-analyze", response_model=QCReport)
async def run_auto_qc(
    request: AutoQCRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Run automatic quality control analysis

    Performs comprehensive QC analysis on a sample and generates detailed report.
    """
    return ai_service.run_auto_qc(request, db)


@router.get("/qc/recommendations/{sample_id}", response_model=List[str])
async def get_qc_recommendations(
    sample_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get QC recommendations for a sample

    Provides actionable recommendations to improve sample quality.
    """
    return ai_service.get_qc_recommendations(sample_id, db)


@router.post("/qc/batch-analyze", response_model=List[QCReport])
async def batch_qc_analysis(
    request: BatchQCRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Batch QC analysis for multiple samples

    Performs QC analysis on multiple samples simultaneously.
    """
    return ai_service.batch_qc_analysis(request, db)


@router.post("/qc/predict-issues", response_model=List[QCIssuePrediction])
async def predict_qc_issues(
    request: PredictQCRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Predict QC issues before sequencing

    Analyzes metadata to predict potential quality issues and prevention strategies.
    """
    return ai_service.predict_qc_issues(request, db)


# ============================================================================
# Anomaly Detection
# ============================================================================

@router.post("/anomaly/detect", response_model=AnomalyDetectionReport)
async def detect_anomalies(
    request: AnomalyDetectionRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Detect anomalies in samples or project

    Uses AI to identify unusual patterns, quality drops, or resource issues.
    """
    return ai_service.detect_anomalies(request, db)


@router.get("/anomaly/{anomaly_id}", response_model=Anomaly)
async def get_anomaly_details(
    anomaly_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get anomaly details

    Retrieves comprehensive information about a specific detected anomaly.
    """
    return ai_service.get_anomaly_details(anomaly_id, db)


@router.post("/anomaly/{anomaly_id}/fix", response_model=AnomalyFixResponse)
async def apply_anomaly_fix(
    anomaly_id: str,
    request: AnomalyFixRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Apply automatic fix for anomaly

    Attempts to automatically resolve detected anomalies when possible.
    """
    return ai_service.apply_anomaly_fix(anomaly_id, request, db)


# ============================================================================
# Smart Sample Grouping
# ============================================================================

@router.post("/grouping/smart-group", response_model=SmartGroupingResult)
async def smart_group_samples(
    request: SmartGroupingRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Perform smart sample grouping

    Uses AI to automatically group samples based on metadata and characteristics.
    """
    return ai_service.smart_group_samples(request, db)


@router.get("/grouping/suggest-comparisons/{project_id}", response_model=List[ComparisonSuggestion])
async def suggest_comparisons(
    project_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Suggest comparison groups

    Recommends meaningful comparisons between sample groups for analysis.
    """
    return ai_service.suggest_comparisons(project_id, db)


@router.post("/grouping/validate", response_model=GroupValidationResult)
async def validate_grouping(
    request: GroupValidationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Validate sample grouping

    Checks sample grouping for statistical validity and potential issues.
    """
    return ai_service.validate_grouping(request, db)


# ============================================================================
# AI Assistant
# ============================================================================

@router.post("/assistant/conversations", response_model=AIAssistantConversation)
async def create_conversation(
    request: CreateConversationRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Create new AI assistant conversation

    Starts a new conversation thread with the AI assistant.
    """
    return ai_service.create_conversation(request, str(current_user.id), db)


@router.get("/assistant/conversations", response_model=List[AIAssistantConversation])
async def list_conversations(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    List all AI assistant conversations

    Retrieves all conversation history for the current user.
    """
    return ai_service.list_conversations(str(current_user.id), db)


@router.get("/assistant/conversations/{conversation_id}", response_model=AIAssistantConversation)
async def get_conversation(
    conversation_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get conversation history

    Retrieves a specific conversation with full message history.
    """
    return ai_service.get_conversation(conversation_id, str(current_user.id), db)


@router.post("/assistant/conversations/{conversation_id}/messages", response_model=AIAssistantMessage)
async def send_assistant_message(
    conversation_id: str,
    request: AssistantMessageRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Send message to AI assistant

    Sends a message and receives AI-generated response.
    """
    return ai_service.send_assistant_message(
        conversation_id,
        request,
        str(current_user.id),
        db
    )


# ============================================================================
# Analysis Insights
# ============================================================================

@router.get("/insights/project/{project_id}", response_model=InsightsReport)
async def get_project_insights(
    project_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get AI-generated insights for a project

    Analyzes project data and provides actionable insights and recommendations.
    """
    return ai_service.get_project_insights(project_id, db)


@router.post("/insights/analyze", response_model=List[AnalysisInsight])
async def get_analysis_insights(
    request: AnalyzeInsightsRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get insights for specific analysis

    Analyzes specific data and generates targeted insights.
    """
    return ai_service.get_analysis_insights(request, db)


@router.put("/insights/{insight_id}/reviewed")
async def mark_insight_reviewed(
    insight_id: str,
    project_id: str = None,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Mark insight as reviewed

    Updates insight status to indicate user has reviewed it.
    """
    ai_service.mark_insight_reviewed(insight_id, project_id, db)
    return {"success": True, "message": "Insight marked as reviewed"}


# ============================================================================
# Predictions
# ============================================================================

@router.post("/predictions/resources", response_model=ResourcePrediction)
async def predict_resources(
    request: ResourcePredictionRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Predict resource requirements for pipeline

    Estimates CPU, memory, storage, and time requirements based on pipeline configuration.
    """
    return ai_service.predict_resources(request, db)


@router.post("/predictions/success", response_model=SuccessPrediction)
async def predict_success(
    analysis_config: Dict[str, Any],
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Predict success probability for analysis

    Analyzes configuration and predicts likelihood of successful completion.
    """
    return ai_service.predict_success(analysis_config, db)


@router.get("/predictions/timeline/{project_id}", response_model=TimelinePrediction)
async def predict_timeline(
    project_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Predict optimal analysis timeline

    Generates timeline predictions with milestones for project completion.
    """
    return ai_service.predict_timeline(project_id, db)


# ============================================================================
# System Status
# ============================================================================

@router.get("/status", response_model=AISystemStatus)
async def get_system_status(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Get AI system status and capabilities

    Returns current operational status and available AI features.
    """
    return ai_service.get_system_status(db)


# ============================================================================
# Feedback & Learning
# ============================================================================

@router.post("/feedback/{recommendation_id}", response_model=FeedbackResponse)
async def submit_feedback(
    recommendation_id: str,
    feedback: FeedbackRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Submit feedback on AI recommendation

    Allows users to provide feedback to improve AI model accuracy.
    """
    return ai_service.submit_feedback(
        recommendation_id,
        feedback,
        str(current_user.id),
        db
    )


@router.post("/feedback/prediction/{prediction_id}", response_model=FeedbackResponse)
async def report_incorrect_prediction(
    prediction_id: str,
    report: IncorrectPredictionRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user)
):
    """
    Report incorrect prediction

    Reports when AI predictions don't match actual outcomes for model improvement.
    """
    return ai_service.report_incorrect_prediction(
        prediction_id,
        report,
        str(current_user.id),
        db
    )
