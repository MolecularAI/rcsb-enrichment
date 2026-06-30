"""Holo structure lookup: resolve InChIKey→CCD, find holo entries, score ligand quality."""

import logging

from .client import RCSBClient, RCSB_DATA, RCSB_SEARCH
from .quality import get_ligand_quality

log = logging.getLogger(__name__)

_MAX_HOLO_ENTRIES = 3  # holo structures evaluated per binder (API cost control)


def _inchikey_to_ccd(client: RCSBClient, inchikey: str) -> str | None:
    """Resolve an InChIKey to a PDB CCD code via the RCSB text_chem search."""
    if not inchikey:
        return None
    payload = {
        "query": {
            "type": "terminal",
            "service": "text_chem",
            "parameters": {
                "attribute": "rcsb_chem_comp_descriptor.InChIKey",
                "operator": "exact_match",
                "value": inchikey,
            },
        },
        "return_type": "mol_definition",
        "request_options": {"paginate": {"start": 0, "rows": 1}},
    }
    data = client.post(RCSB_SEARCH, payload)
    results = (data or {}).get("result_set") or []
    return results[0]["identifier"] if results else None


def _find_holo_entries(
    client: RCSBClient,
    uniprot_id: str,
    ccd_code: str,
    max_rows: int = _MAX_HOLO_ENTRIES,
) -> list:
    """Find PDB entries that contain both a given UniProt target and a co-crystallised ligand."""
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                                     ".reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_nonpolymer_entity_container_identifiers"
                                     ".nonpolymer_comp_id",
                        "operator": "exact_match",
                        "value": ccd_code,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows}},
    }
    data = client.post(RCSB_SEARCH, payload)
    return [r["identifier"] for r in (data or {}).get("result_set") or []]


def get_holo_ligand_quality(
    client: RCSBClient,
    binders: list,
    uniprot_ids: list,
) -> list:
    """For each known binder, find holo PDB entries and score ligand binding quality.

    Steps per binder:
    1. Resolve InChIKey → CCD code (skip if unresolvable)
    2. Search for holo entries: UniProt + CCD code
    3. Run get_ligand_quality() on each holo entry (up to _MAX_HOLO_ENTRIES)
    4. Attach results to the binder dict

    Returns a list of dicts, one per binder that has at least one holo entry found.
    """
    if not uniprot_ids:
        return []

    results = []
    uniprot_id = uniprot_ids[0]

    for b in binders:
        ccd = b.get("chem_comp_id")
        if not ccd:
            ccd = _inchikey_to_ccd(client, b.get("inchikey") or "")
        if not ccd:
            log.debug("No CCD resolvable for binder %s — skipping holo search", b.get("name"))
            continue

        holo_ids = _find_holo_entries(client, uniprot_id, ccd)
        if not holo_ids:
            log.info("No holo entries found for UniProt %s + CCD %s", uniprot_id, ccd)
            continue

        log.info("Binder %s (CCD %s): %d holo entr%s found: %s",
                 b.get("name") or ccd, ccd, len(holo_ids),
                 "y" if len(holo_ids) == 1 else "ies", ",".join(holo_ids))

        holo_metrics = []
        for holo_pdb in holo_ids:
            entry = client.get(f"{RCSB_DATA}/entry/{holo_pdb}")
            if not entry:
                continue
            np_ids = (
                (entry.get("rcsb_entry_container_identifiers") or {})
                .get("non_polymer_entity_ids") or []
            )
            metrics = get_ligand_quality(client, holo_pdb, np_ids)
            holo_metrics.extend(
                {**m, "holo_pdb": holo_pdb}
                for m in metrics if m.get("ligand_id") == ccd
            )

        if holo_metrics:
            results.append({
                "binder_name": b.get("name"),
                "binder_inchikey": b.get("inchikey"),
                "ccd_code": ccd,
                "holo_pdb_ids": ",".join(holo_ids),
                "holo_ligand_metrics": holo_metrics,
            })

    return results
