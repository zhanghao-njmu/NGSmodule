"""
Base AI provider interface.

All concrete providers (Anthropic, OpenAI, local mock) implement this contract
so the rest of the system can swap them at runtime.
"""

from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional


class AIProviderError(Exception):
    """Raised when an AI provider fails to produce a response."""


class AIProvider(ABC):
    """Abstract AI provider.

    Concrete providers must implement at least `chat()` for free-form
    conversational responses. The other helpers have sensible defaults
    that build on top of `chat()`.
    """

    name: str = "base"

    @abstractmethod
    def chat(
        self,
        messages: List[Dict[str, str]],
        system: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 1024,
        **kwargs: Any,
    ) -> str:
        """Send a list of messages and return a single text response.

        Each message is a dict with keys "role" (user|assistant|system)
        and "content" (str).
        """

    def complete(self, prompt: str, **kwargs: Any) -> str:
        """Single-prompt convenience wrapper."""
        return self.chat([{"role": "user", "content": prompt}], **kwargs)

    def is_available(self) -> bool:
        """Whether the provider is configured and reachable."""
        return True

    def info(self) -> Dict[str, Any]:
        """Return human-readable provider info."""
        return {"name": self.name, "available": self.is_available()}
