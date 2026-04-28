"""
Tests for the data-downloads + vendor-credentials API surface.

Vendor adapters are stubbed via monkeypatch so tests don't shell out to
the real lcbio Java daemon. We focus on the contract: routing, schema,
DB persistence, encryption.
"""
from pathlib import Path

import pytest

from app.core.security import encrypt_secret, decrypt_secret, try_decrypt_secret
from app.services.data_provider.base import (
    DataProviderBase,
    SessionInfo,
    Progress,
)


# ---------------- crypto helpers ----------------

def test_encrypt_decrypt_roundtrip():
    """Plain → encrypt → decrypt should give back the same string."""
    plaintext = "p@ssw0rd! 中文测试"
    cipher = encrypt_secret(plaintext)
    assert cipher != plaintext
    assert decrypt_secret(cipher) == plaintext


def test_try_decrypt_returns_none_for_garbage():
    """try_decrypt_secret should not raise on a malformed token."""
    assert try_decrypt_secret("not-a-fernet-token") is None


def test_two_encryptions_give_different_ciphers():
    """Fernet uses a random IV — same plaintext, different ciphertext."""
    a = encrypt_secret("same")
    b = encrypt_secret("same")
    assert a != b
    assert decrypt_secret(a) == decrypt_secret(b) == "same"


# ---------------- vendor credential CRUD ----------------

def test_credential_create_then_list_masked(client, auth_headers):
    """Creating a credential returns masked email; password never echoed."""
    r = client.post(
        "/api/v1/vendor-credentials",
        headers=auth_headers,
        json={
            "vendor": "lc_bio",
            "label": "default",
            "email": "alice@example.com",
            "password": "topsecret",
        },
    )
    assert r.status_code == 201, r.text
    body = r.json()
    assert body["vendor"] == "lc_bio"
    assert body["email_preview"].startswith("al")
    assert "@example.com" in body["email_preview"]
    assert "topsecret" not in r.text  # raw password must not leak

    r = client.get("/api/v1/vendor-credentials", headers=auth_headers)
    assert r.status_code == 200
    items = r.json()["items"]
    assert len(items) == 1
    assert items[0]["vendor"] == "lc_bio"


def test_credential_overwrite_on_same_vendor_label(client, auth_headers):
    """Re-POSTing the same (vendor, label) updates instead of creating."""
    base = {"vendor": "lc_bio", "label": "shared", "email": "old@example.com", "password": "old-pw"}
    r1 = client.post("/api/v1/vendor-credentials", headers=auth_headers, json=base)
    assert r1.status_code == 201
    cid = r1.json()["id"]

    base["email"] = "new@example.com"
    base["password"] = "new-pw"
    r2 = client.post("/api/v1/vendor-credentials", headers=auth_headers, json=base)
    assert r2.status_code == 201
    assert r2.json()["id"] == cid  # same row, updated in place

    r = client.get("/api/v1/vendor-credentials", headers=auth_headers)
    items = r.json()["items"]
    assert len(items) == 1


def test_credential_delete(client, auth_headers):
    r = client.post(
        "/api/v1/vendor-credentials",
        headers=auth_headers,
        json={"vendor": "lc_bio", "email": "to-delete@example.com", "password": "x"},
    )
    cid = r.json()["id"]

    r2 = client.delete(f"/api/v1/vendor-credentials/{cid}", headers=auth_headers)
    assert r2.status_code == 200

    r3 = client.get("/api/v1/vendor-credentials", headers=auth_headers)
    assert r3.json()["total"] == 0


# ---------------- vendor list ----------------

def test_supported_vendors(client, auth_headers):
    """Both lc_bio (real) and novogene (stub) should advertise."""
    r = client.get("/api/v1/data-downloads/vendors", headers=auth_headers)
    assert r.status_code == 200
    vendors = r.json()
    assert "lc_bio" in vendors
    assert "novogene" in vendors


# ---------------- download job lifecycle (stubbed provider) ----------------

class _StubProvider(DataProviderBase):
    """In-memory provider for tests — no Java daemon required."""
    name = "lc_bio"
    _session = False
    _log_path: Path

    def session_status(self):
        return SessionInfo(active=self._session, pid=99 if self._session else None)

    def login(self, _email, _password):
        self._session = True
        return self.session_status()

    def logout(self):
        self._session = False

    def start_download(self, source_path, dest_dir):
        dest_dir = Path(dest_dir)
        dest_dir.mkdir(parents=True, exist_ok=True)
        log = dest_dir / f"{Path(source_path).name}_下载日志_test.txt"
        log.write_text("已开始下载，请耐心等待......\n")
        self._log_path = log
        return log

    def parse_progress(self, log_path):
        text = Path(log_path).read_text()
        if "下载已完成" in text:
            return Progress(status="completed", percent=100.0)
        if "已下载：" in text:
            # Stub returns 50% to verify the field flows through
            return Progress(status="running", percent=50.0, bytes_downloaded=512 * 1024 * 1024)
        return Progress(status="running", percent=0.0)


@pytest.fixture
def stubbed_provider(monkeypatch):
    """Replace lc_bio provider in factory with the in-memory stub."""
    stub = _StubProvider()
    # Patch every callsite (factory + the modules that imported it).
    monkeypatch.setattr("app.services.data_provider.factory.get_provider", lambda _v: stub)
    monkeypatch.setattr("app.services.data_download_service.get_provider", lambda _v: stub)
    monkeypatch.setattr("app.services.data_provider.get_provider", lambda _v: stub, raising=False)
    return stub


def test_download_job_create_and_list(client, auth_headers, stubbed_provider, tmp_path):
    """Open session → create job → list returns the job."""
    r = client.post(
        "/api/v1/data-downloads/sessions",
        headers=auth_headers,
        json={"vendor": "lc_bio", "email": "x@y", "password": "z"},
    )
    assert r.status_code == 201, r.text
    assert r.json()["active"] is True

    r = client.post(
        "/api/v1/data-downloads/jobs",
        headers=auth_headers,
        json={
            "vendor": "lc_bio",
            "source_path": "PROJ-1/2026/Data.tar",
            "dest_path": str(tmp_path),
        },
    )
    assert r.status_code == 201, r.text
    job = r.json()
    assert job["status"] == "running"
    assert job["vendor"] == "lc_bio"

    # Simulate the download making progress
    Path(stubbed_provider._log_path).write_text(
        "已开始下载，请耐心等待......\n"
        "PROJ-1/2026/Data.tar ：已下载：512.00 MB  50.00%      2026-04-28-10:00:00.000\n"
    )

    r = client.get(f"/api/v1/data-downloads/jobs/{job['id']}", headers=auth_headers)
    assert r.status_code == 200
    assert r.json()["progress_pct"] == 50.0


def test_download_job_completion(client, auth_headers, stubbed_provider, tmp_path):
    """Completion marker flips status to 'completed' and freezes progress at 100."""
    client.post(
        "/api/v1/data-downloads/sessions",
        headers=auth_headers,
        json={"vendor": "lc_bio", "email": "x@y", "password": "z"},
    )
    r = client.post(
        "/api/v1/data-downloads/jobs",
        headers=auth_headers,
        json={
            "vendor": "lc_bio",
            "source_path": "PROJ-1/2026/Data.tar",
            "dest_path": str(tmp_path),
        },
    )
    job = r.json()

    Path(stubbed_provider._log_path).write_text(
        "PROJ-1/2026/Data.tar ： 下载已完成！      2026-04-28-10:00:00.000\n"
    )
    r = client.get(f"/api/v1/data-downloads/jobs/{job['id']}", headers=auth_headers)
    assert r.json()["status"] == "completed"
    assert r.json()["progress_pct"] == 100.0


def test_session_login_validation_error(client, auth_headers):
    """POST /sessions with neither credential_id nor email/password rejects."""
    r = client.post(
        "/api/v1/data-downloads/sessions",
        headers=auth_headers,
        json={"vendor": "lc_bio"},
    )
    assert r.status_code == 400
    assert "credential_id" in r.json()["detail"]


# ---------------- lcbio log parser ----------------

def test_lcbio_parse_progress_line():
    """Sanity-check the regex against actual log samples."""
    from app.services.data_provider.lc_bio import LcBioProvider

    log = Path(__file__).parent / "_lcbio_sample.txt"
    log.write_text(
        "已开始下载，请耐心等待......\n"
        "X ：下载中......       2026-04-27-16:20:54.820\n"
        "X ：已下载：1.03 GB  8.15%      2026-04-27-16:22:32.392\n"
    )
    p = LcBioProvider().parse_progress(log)
    assert p.status == "running"
    assert p.percent == pytest.approx(8.15)
    assert p.bytes_downloaded is not None and p.bytes_downloaded > 1024 ** 3 * 0.99


def test_lcbio_parse_progress_completed(tmp_path):
    from app.services.data_provider.lc_bio import LcBioProvider

    log = tmp_path / "done.txt"
    log.write_text("X ： 下载已完成！      2026-04-27-16:38:47.328\n")
    p = LcBioProvider().parse_progress(log)
    assert p.status == "completed"
    assert p.percent == 100.0
