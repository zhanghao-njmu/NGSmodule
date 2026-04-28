"""
Mock AI provider used in development and tests.

Returns plausible canned responses without contacting any external service.
"""

import random
from typing import Any, Dict, List, Optional

from app.services.ai_providers.base import AIProvider

_RESPONSES = [
    "Based on your data, I recommend using a quality threshold of 20 for optimal results.",
    "I've analyzed your pipeline configuration and found potential optimization opportunities.",
    "The anomaly you're experiencing is likely due to adapter contamination. Try running Cutadapt.",
    "Your samples show good clustering by condition. Proceed with differential expression analysis.",
    "I suggest increasing the thread count to 16 for better performance on your dataset.",
]


class MockAIProvider(AIProvider):
    name = "mock"

    def chat(
        self,
        messages: List[Dict[str, str]],
        system: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 1024,
        **kwargs: Any,
    ) -> str:
        # Pick a deterministic-ish response keyed off the last user message
        if messages:
            last = messages[-1].get("content", "")
            idx = abs(hash(last)) % len(_RESPONSES)
            return _RESPONSES[idx]
        return random.choice(_RESPONSES)
