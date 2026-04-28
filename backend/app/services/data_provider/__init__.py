"""
Data provider adapters for sequencing vendors.

Each vendor (联川 / 诺禾致源 / ...) implements `DataProviderBase` to
expose a uniform login + download interface to the rest of the app.
"""
from app.services.data_provider.base import (
    DataProviderBase,
    SessionExistsError,
    NoSessionError,
    LoginFailedError,
    DownloadFailedError,
)
from app.services.data_provider.factory import get_provider

__all__ = [
    "DataProviderBase",
    "SessionExistsError",
    "NoSessionError",
    "LoginFailedError",
    "DownloadFailedError",
    "get_provider",
]
