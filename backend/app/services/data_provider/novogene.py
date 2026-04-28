"""
Novogene 诺禾致源 data download adapter — STUB.

The actual integration depends on which delivery method Novogene provides
for a given customer (lcbio-style daemon? sftp? web download token? OBS
presigned URL?). This stub keeps the registry slot reserved so the
backend exposes the vendor in `/data-downloads/vendors` and the frontend
can render a disabled-with-tooltip option.

To finish:
  1. Identify the actual delivery channel from Novogene's docs/customer
     account (e.g. ftp.novogene.com / s3 bucket / their own CLI).
  2. Replace each method below with the corresponding logic, mirroring
     `lc_bio.py`.
  3. Register in `factory.py` once tested end-to-end.
"""
from pathlib import Path

from app.services.data_provider.base import (
    DataProviderBase,
    SessionInfo,
    Progress,
    LoginFailedError,
)


class NovogeneProvider(DataProviderBase):
    name = "novogene"

    def session_status(self) -> SessionInfo:
        return SessionInfo(active=False)

    def login(self, email: str, password: str) -> SessionInfo:
        raise LoginFailedError(
            "Novogene adapter not implemented yet. "
            "See backend/app/services/data_provider/novogene.py for the integration TODO."
        )

    def logout(self) -> None:
        return

    def start_download(self, source_path: str, dest_dir: Path) -> Path:
        raise NotImplementedError("Novogene adapter not implemented yet.")

    def parse_progress(self, log_path: Path) -> Progress:
        return Progress(status="failed", percent=0.0, error_message="adapter not implemented")
