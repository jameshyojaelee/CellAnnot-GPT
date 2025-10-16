from __future__ import annotations

import sys
from types import SimpleNamespace

pytest_plugins = ("pytest_asyncio", "pytest_asyncio.plugin")


class _DummyOpenAI:
    def __init__(self, *_, **__) -> None:
        self.chat = SimpleNamespace(completions=SimpleNamespace(create=lambda **kwargs: None))


sys.modules.setdefault("openai", SimpleNamespace(OpenAI=_DummyOpenAI))
