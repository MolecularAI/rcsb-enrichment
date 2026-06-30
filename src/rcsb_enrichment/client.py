"""HTTP client and TLS setup for RCSB/PDBe API requests."""

import logging
import os
import re
import ssl
import subprocess
import tempfile
import time
from typing import Any

log = logging.getLogger(__name__)

import certifi
import requests

RCSB_DATA = "https://data.rcsb.org/rest/v1/core"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBE_API = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites"

_PEM_RE = re.compile(
    rb"-----BEGIN CERTIFICATE-----.*?-----END CERTIFICATE-----", re.DOTALL
)


def _cert_accepted_by_ssl(pem_block: bytes) -> bool:
    """Return True if OpenSSL accepts this cert as a CA in a verify context."""
    ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        ctx.load_verify_locations(cadata=pem_block.decode("ascii", errors="replace"))
        return True
    except ssl.SSLError:
        return False


def _build_ca_bundle() -> str:
    """Combine macOS system keychain certs with certifi's bundle.

    On networks with TLS-inspecting proxies the system keychain contains
    the proxy's root CA which certifi's bundle does not include.
    Each keychain cert is validated against the running OpenSSL before
    inclusion so that non-compliant certs (e.g. Basic Constraints not
    marked critical, rejected by Python 3.14+) are silently dropped.
    Returns a path to a combined PEM file (written once per process).
    """
    keychains = [
        "/Library/Keychains/System.keychain",
        "/System/Library/Keychains/SystemRootCertificates.keychain",
    ]
    good_pem_blocks: list[bytes] = []
    for kc in keychains:
        if os.path.exists(kc):
            r = subprocess.run(
                ["security", "find-certificate", "-a", "-p", kc],
                capture_output=True,
            )
            if r.returncode == 0:
                for block in _PEM_RE.findall(r.stdout):
                    if _cert_accepted_by_ssl(block):
                        good_pem_blocks.append(block)

    with open(certifi.where(), "rb") as f:
        good_pem_blocks.extend(_PEM_RE.findall(f.read()))

    tmp = tempfile.NamedTemporaryFile(mode="wb", suffix=".pem", delete=False, prefix="rcsb_ca_")
    tmp.write(b"\n".join(good_pem_blocks))
    tmp.close()
    return tmp.name


_CA_BUNDLE = _build_ca_bundle()


class RCSBClient:
    def __init__(self, delay: float = 0.1):
        self._delay = delay
        self._session = requests.Session()
        self._session.headers.update({"Accept": "application/json"})
        self._session.verify = _CA_BUNDLE

    def get(self, url: str) -> Any:
        time.sleep(self._delay)
        r = self._session.get(url, timeout=30)
        if r.status_code in (404, 204):
            return None
        if r.status_code == 429:
            wait = int(r.headers.get("Retry-After", 5))
            log.warning("Rate limited; waiting %ds", wait)
            time.sleep(wait)
            return self.get(url)
        r.raise_for_status()
        if not r.content or not r.text.strip():
            return None
        return r.json()

    def post(self, url: str, payload: dict) -> Any:
        time.sleep(self._delay)
        r = self._session.post(url, json=payload, timeout=30)
        if r.status_code in (404, 204):
            return None
        if r.status_code == 429:
            wait = int(r.headers.get("Retry-After", 5))
            log.warning("Rate limited; waiting %ds", wait)
            time.sleep(wait)
            return self.post(url, payload)
        r.raise_for_status()
        if not r.content or not r.text.strip():
            return None
        return r.json()
