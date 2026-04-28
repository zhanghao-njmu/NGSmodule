"""
Anthropic AI provider.

Uses the Anthropic Python SDK with prompt caching enabled. Requires
`pip install anthropic` and ANTHROPIC_API_KEY in settings.
"""
import logging
from typing import Any, Dict, List, Optional

from app.core.config import settings
from app.services.ai_providers.base import AIProvider, AIProviderError

logger = logging.getLogger(__name__)


class AnthropicProvider(AIProvider):
    name = "anthropic"

    # Default model — operators can override via settings.ANTHROPIC_MODEL.
    # Note: this string is the literal API parameter; do not paraphrase.
    default_model: str = "claude-sonnet-4-6"

    def __init__(self):
        self._client = None
        self._api_key = getattr(settings, "ANTHROPIC_API_KEY", None)
        self._model = getattr(settings, "ANTHROPIC_MODEL", self.default_model)

    def _get_client(self):
        if self._client is None:
            if not self._api_key:
                raise AIProviderError(
                    "Anthropic provider selected but ANTHROPIC_API_KEY is not configured"
                )
            try:
                from anthropic import Anthropic
            except ImportError as exc:
                raise AIProviderError(
                    "anthropic SDK is not installed. Run: pip install anthropic"
                ) from exc
            self._client = Anthropic(api_key=self._api_key)
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

        # Anthropic API uses a separate `system` parameter and rejects
        # role="system" inside `messages`. Split if needed.
        sys_prompt = system
        anthropic_messages = []
        for m in messages:
            role = m.get("role", "user")
            content = m.get("content", "")
            if role == "system" and not sys_prompt:
                sys_prompt = content
            elif role in ("user", "assistant"):
                anthropic_messages.append({"role": role, "content": content})

        try:
            kwargs_api = {
                "model": self._model,
                "max_tokens": max_tokens,
                "temperature": temperature,
                "messages": anthropic_messages,
            }
            if sys_prompt:
                # Prompt caching: cache the system prompt so subsequent
                # requests in the same conversation re-use the parsed tokens
                kwargs_api["system"] = [
                    {
                        "type": "text",
                        "text": sys_prompt,
                        "cache_control": {"type": "ephemeral"},
                    }
                ]

            response = client.messages.create(**kwargs_api)
            blocks = response.content or []
            return "".join(getattr(b, "text", "") for b in blocks)
        except Exception as exc:
            logger.exception("Anthropic API call failed")
            raise AIProviderError(f"Anthropic API error: {exc}") from exc
