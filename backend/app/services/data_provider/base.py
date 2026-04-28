"""
Abstract base class for vendor data-download providers.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# ---------------- exceptions ----------------


class DataProviderError(Exception):
    """Base for provider errors. Surfaced as 4xx/5xx by the API layer."""


class SessionExistsError(DataProviderError):
    """Vendor daemon is already running, cannot start a new session."""


class NoSessionError(DataProviderError):
    """A vendor session is required but none is active."""


class LoginFailedError(DataProviderError):
    """Credentials rejected by the vendor."""


class DownloadFailedError(DataProviderError):
    """Vendor reported the download failed (bad path, network, ...)."""


# ---------------- value objects ----------------


@dataclass
class SessionInfo:
    active: bool
    pid: Optional[int] = None


@dataclass
class Progress:
    """Snapshot of an in-flight download parsed from the vendor's log file."""

    status: str  # "pending" | "running" | "completed" | "failed"
    percent: float  # 0.0 – 100.0
    bytes_downloaded: Optional[int] = None
    file_size: Optional[int] = None
    error_message: Optional[str] = None


# ---------------- provider base ----------------


class DataProviderBase(ABC):
    """Vendor-agnostic contract.

    Concrete adapters wrap whatever native CLI / SDK the vendor ships.
    The session is process-level (one daemon per vendor); concurrent
    downloads inside one session are vendor-dependent.
    """

    name: str  # short id, e.g. "lc_bio"

    @abstractmethod
    def session_status(self) -> SessionInfo:
        """Whether the vendor daemon is currently running."""

    @abstractmethod
    def login(self, email: str, password: str) -> SessionInfo:
        """Start the vendor daemon and authenticate.

        Raises:
            SessionExistsError: daemon already running.
            LoginFailedError: wrong credentials.
        """

    @abstractmethod
    def logout(self) -> None:
        """Tear down the vendor daemon (kills all in-flight downloads)."""

    @abstractmethod
    def start_download(self, source_path: str, dest_dir: Path) -> Path:
        """Trigger a download. Returns the absolute path of the vendor log
        file that will be tailed for progress events.

        Raises:
            NoSessionError: daemon not running.
            DownloadFailedError: vendor rejected the request synchronously.
        """

    @abstractmethod
    def parse_progress(self, log_path: Path) -> Progress:
        """Read the vendor log file and return the latest progress."""
