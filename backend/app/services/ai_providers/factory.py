"""
AI provider factory.

Selects the configured provider based on settings.AI_PROVIDER, falling back
to the mock provider when the requested backend is unavailable.
"""
import logging
from typing import Optional

from app.core.config import settings
from app.services.ai_providers.base import AIProvider, AIProviderError
from app.services.ai_providers.mock import MockAIProvider

logger = logging.getLogger(__name__)

_cached: Optional[AIProvider] = None


def _build_provider(name: str) -> AIProvider:
    name = (name or "mock").lower()

    if name == "claude":
        from app.services.ai_providers.claude import ClaudeProvider
        return ClaudeProvider()

    if name == "openai":
        from app.services.ai_providers.openai_provider import OpenAIProvider
        return OpenAIProvider()

    if name == "mock":
        return MockAIProvider()

    raise AIProviderError(f"Unknown AI provider: {name}")


def get_ai_provider(refresh: bool = False) -> AIProvider:
    """Return the configured AI provider, falling back to mock on errors."""
    global _cached
    if _cached is not None and not refresh:
        return _cached

    name = getattr(settings, "AI_PROVIDER", "mock")
    try:
        provider = _build_provider(name)
        if not provider.is_available():
            logger.warning(
                f"AI provider '{name}' is configured but not available; "
                "falling back to mock provider"
            )
            provider = MockAIProvider()
    except AIProviderError as exc:
        logger.warning(f"AI provider initialization failed ({exc}); using mock")
        provider = MockAIProvider()

    _cached = provider
    logger.info(f"AI provider in use: {provider.info()}")
    return provider
