"""
Vendor → provider lookup.
"""
from typing import Dict, Type

from app.services.data_provider.base import DataProviderBase
from app.services.data_provider.lc_bio import LcBioProvider


_REGISTRY: Dict[str, Type[DataProviderBase]] = {
    "lc_bio": LcBioProvider,
    # "novogene": NovogeneProvider,   # Phase 6
}


def get_provider(vendor: str) -> DataProviderBase:
    """Return a provider instance for the given vendor id.

    Each call returns a fresh instance (cheap), letting callers avoid
    sharing process state across requests.
    """
    cls = _REGISTRY.get(vendor)
    if cls is None:
        raise ValueError(
            f"Unknown vendor {vendor!r}. Supported: {list(_REGISTRY)}"
        )
    return cls()


def list_vendors() -> list[str]:
    return list(_REGISTRY.keys())
