"""
AI provider abstractions.

Allows swapping the AI backend (mock, Anthropic, OpenAI, local model) without
changing the consumer code in `ai_service.py`.
"""
from app.services.ai_providers.base import AIProvider, AIProviderError
from app.services.ai_providers.factory import get_ai_provider

__all__ = ["AIProvider", "AIProviderError", "get_ai_provider"]
