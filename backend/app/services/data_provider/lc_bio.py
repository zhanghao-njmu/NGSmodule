"""
联川生物云平台 (lcbio) data download adapter.

Wraps the bash CLI that ships with their `lcbio_download.run` installer:
  - bin/run <email> <password>   — start the Java daemon and login
  - bin/download <obs_path> <dest_dir>  — trigger an async download
  - bin/stop                     — kill the daemon

The Java daemon (`obs_service-...-SNAPSHOT.jar`) is process-level singleton;
this adapter expects it to be reachable via the shared `LC_BIO_BIN_DIR`
mountpoint when the backend runs containerised.
"""
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Optional

from app.core.config import settings
from app.services.data_provider.base import (
    DataProviderBase,
    SessionInfo,
    Progress,
    SessionExistsError,
    NoSessionError,
    LoginFailedError,
    DownloadFailedError,
)


# JAR name used for `ps` detection
JAR_MARKER = "obs_service-0.0.1-SNAPSHOT.jar"

# Progress line examples (Chinese):
#   <obs> ：已下载：50.00 MB  0.39%      2026-04-27-16:20:56.441
#   <obs> ：已下载：1.03 GB  8.15%       2026-04-27-16:22:32.392
_PROGRESS_RE = re.compile(
    r"已下载：([\d.]+)\s*(MB|GB|TB|KB)\s+([\d.]+)%"
)
_DONE_MARKER = "下载已完成"
_FAIL_MARKERS = ("下载失败", "Error", "失败", "网络异常")
_RUNNING_MARKER = "下载中"


def _to_bytes(amount: float, unit: str) -> int:
    multipliers = {"KB": 1024, "MB": 1024 ** 2, "GB": 1024 ** 3, "TB": 1024 ** 4}
    return int(amount * multipliers.get(unit, 1))


class LcBioProvider(DataProviderBase):
    name = "lc_bio"

    def __init__(self, bin_dir: Optional[Path] = None):
        # Resolve at construction so tests can override; settings provides default.
        self.bin_dir = Path(bin_dir or getattr(settings, "LC_BIO_BIN_DIR", "/home/zhanghao/programs/lcbio/linux_download/bin"))
        self.run_bin = self.bin_dir / "run"
        self.download_bin = self.bin_dir / "download"
        self.stop_bin = self.bin_dir / "stop"

    # ------------ session ------------

    def session_status(self) -> SessionInfo:
        pid = self._daemon_pid()
        return SessionInfo(active=pid is not None, pid=pid)

    def login(self, email: str, password: str) -> SessionInfo:
        if self._daemon_pid() is not None:
            raise SessionExistsError("lcbio daemon already running; logout first")

        # `run` is interactive — it greps the log for ResponseContent and then
        # waits for an Enter keypress. We close stdin so the read returns EOF
        # and the script exits, then verify the daemon spawned in the
        # background.
        proc = subprocess.run(
            [str(self.run_bin), email, password],
            capture_output=True,
            text=True,
            timeout=60,
            stdin=subprocess.DEVNULL,
        )
        # Wait briefly for daemon to come up (the run script backgrounds it).
        for _ in range(30):
            if self._daemon_pid() is not None:
                break
            time.sleep(0.2)
        else:
            stderr = (proc.stderr or proc.stdout or "").strip()
            raise LoginFailedError(f"lcbio daemon did not start: {stderr[:300]}")

        # The Java service truncates linuxLogin.txt to 0 bytes once it has read
        # the credentials; if it's still non-empty after a few seconds we
        # likely have a credential problem.
        login_file = self.bin_dir.parent / "lib" / ".obs_service" / "linuxLogin.txt"
        for _ in range(50):
            try:
                if login_file.stat().st_size == 0:
                    break
            except FileNotFoundError:
                break
            time.sleep(0.1)

        return self.session_status()

    def logout(self) -> None:
        if self._daemon_pid() is None:
            return
        subprocess.run(
            [str(self.stop_bin)],
            capture_output=True,
            text=True,
            timeout=15,
            stdin=subprocess.DEVNULL,
        )

    # ------------ download ------------

    def start_download(self, source_path: str, dest_dir: Path) -> Path:
        if self._daemon_pid() is None:
            raise NoSessionError("lcbio daemon not running; login first")

        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)

        # `download` is also interactive; same DEVNULL trick + short timeout
        # because the actual download runs inside the Java daemon.
        try:
            subprocess.run(
                [str(self.download_bin), source_path, str(dest_dir)],
                capture_output=True,
                text=True,
                timeout=15,
                stdin=subprocess.DEVNULL,
            )
        except subprocess.TimeoutExpired:
            # Expected — printlog.sh hangs on `read -p`. The download has
            # already been queued in linuxDown.txt by then.
            pass
        except subprocess.CalledProcessError as exc:
            raise DownloadFailedError(f"lcbio download command failed: {exc.stderr or exc.stdout}")

        # Vendor writes one log per download named:
        #   <basename>_下载日志_YYYY-MM-DD-HH-MM-SS.txt
        basename = source_path.rstrip("/").split("/")[-1]
        log = self._wait_for_log(dest_dir, basename, timeout=10)
        if log is None:
            raise DownloadFailedError(
                f"lcbio download log not produced under {dest_dir} for {basename!r}; "
                "check daemon health or source path"
            )
        return log

    def parse_progress(self, log_path: Path) -> Progress:
        if not log_path.exists():
            return Progress(status="pending", percent=0.0)

        try:
            # Tail the last ~50 lines; progress lines are short and the file
            # can grow into MB-scale, so avoid reading the whole thing.
            with open(log_path, "rb") as f:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                f.seek(max(0, size - 8192))
                tail = f.read().decode("utf-8", errors="replace")
        except OSError:
            return Progress(status="pending", percent=0.0)

        lines = [ln for ln in tail.splitlines() if ln.strip()]
        if not lines:
            return Progress(status="pending", percent=0.0)

        last = lines[-1]
        if _DONE_MARKER in last:
            return Progress(status="completed", percent=100.0)
        if any(m in last for m in _FAIL_MARKERS):
            return Progress(status="failed", percent=0.0, error_message=last.strip())

        m = _PROGRESS_RE.search(last)
        if m:
            amount, unit, pct = float(m.group(1)), m.group(2), float(m.group(3))
            return Progress(
                status="running",
                percent=pct,
                bytes_downloaded=_to_bytes(amount, unit),
            )

        if _RUNNING_MARKER in last:
            return Progress(status="running", percent=0.0)

        return Progress(status="running", percent=0.0)

    # ------------ helpers ------------

    @staticmethod
    def _daemon_pid() -> Optional[int]:
        try:
            out = subprocess.check_output(["ps", "-eo", "pid,command"], text=True, timeout=5)
        except (subprocess.SubprocessError, FileNotFoundError):
            return None
        for line in out.splitlines():
            if JAR_MARKER in line:
                try:
                    return int(line.strip().split(None, 1)[0])
                except (ValueError, IndexError):
                    continue
        return None

    @staticmethod
    def _wait_for_log(dest_dir: Path, basename: str, timeout: float) -> Optional[Path]:
        deadline = time.time() + timeout
        pattern = f"{basename}_下载日志_*.txt"
        # Log filename may be 'Data.tar_下载日志_*.txt' or just '下载日志_*.txt'
        # depending on lcbio version — try both.
        fallback = "下载日志_*.txt"
        while time.time() < deadline:
            matches = sorted(dest_dir.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
            if not matches:
                matches = sorted(dest_dir.glob(fallback), key=lambda p: p.stat().st_mtime, reverse=True)
            if matches:
                return matches[0]
            time.sleep(0.5)
        return None
