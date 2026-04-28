"""
OpenAI provider.

Requires `pip install openai` and OPENAI_API_KEY in settings.
"""

import logging
from typing import Any, Dict, List, Optional

from app.core.config import settings
from app.services.ai_providers.base import AIProvider, AIProviderError

logger = logging.getLogger(__name__)


class OpenAIProvider(AIProvider):
    name = "openai"
    default_model: str = "gpt-4o-mini"

    def __init__(self):
        self._client = None
        self._api_key = getattr(settings, "OPENAI_API_KEY", None)
        self._model = getattr(settings, "OPENAI_MODEL", self.default_model)

    def _get_client(self):
        if self._client is None:
            if not self._api_key:
                raise AIProviderError("OpenAI provider selected but OPENAI_API_KEY is not configured")
            try:
                from openai import OpenAI
            except ImportError as exc:
                raise AIProviderError("openai SDK is not installed. Run: pip install openai") from exc
            self._client = OpenAI(api_key=self._api_key)
        return self._client

    def is_available(self) -> bool:
        return bool(self._api_key)

    def info(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "model": self._model,
            "available": self.is_available(),
        }

    def chat(
        self,
        messages: List[Dict[str, str]],
        system: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 1024,
        **kwargs: Any,
    ) -> str:
        client = self._get_client()

        api_messages = list(messages)
        if system and not any(m.get("role") == "system" for m in api_messages):
            api_messages = [{"role": "system", "content": system}] + api_messages

        try:
            response = client.chat.completions.create(
                model=self._model,
                messages=api_messages,
                temperature=temperature,
                max_tokens=max_tokens,
            )
            return response.choices[0].message.content or ""
        except Exception as exc:
            logger.exception("OpenAI API call failed")
            raise AIProviderError(f"OpenAI API error: {exc}") from exc
