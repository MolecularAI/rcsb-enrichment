"""PDBe binding site data (SIFTS-curated)."""

from .client import RCSBClient, PDBE_API


def get_pdbe_binding_sites(client: RCSBClient, pdb_id: str) -> list:
    data = client.get(f"{PDBE_API}/{pdb_id.lower()}")
    if not data:
        return []

    entry_data = data.get(pdb_id.lower()) or []
    return [
        {
            "description": site.get("site_description") or site.get("details", ""),
            "source": "PDBe-KB/SIFTS",
        }
        for site in entry_data
    ]
