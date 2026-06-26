"""
RCSB PDB enrichment script.

Reads PDB IDs from a CSV, queries RCSB REST APIs, and writes an enriched CSV
with quality metrics, related PDB entries, and binding site annotations.

Usage:
    python rcsb_enrichment.py --input proteins.csv --output enriched.csv
    python rcsb_enrichment.py --input proteins.csv --output enriched.csv \
        --pdb-col pdb_id --uniprot-col uniprot_id \
        --seq-identity 0.9 --max-related 25 --delay 0.1
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import tempfile
import time
from typing import Any

import certifi
import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

RCSB_DATA = "https://data.rcsb.org/rest/v1/core"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBE_API = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites"

CONTACT_DIST_CUTOFF = 4.0  # Å — ligand–residue contact threshold

_PDB_COL_ALIASES = {"pdb_id", "pdbid", "pdb", "pdb id"}


def _build_ca_bundle() -> str:
    """Combine macOS system keychain certs with certifi's bundle.

    On networks with TLS-inspecting proxies the system keychain contains
    the proxy's root CA which certifi's bundle does not include.
    Returns a path to a combined PEM file (written once per process).
    """
    keychains = [
        "/Library/Keychains/System.keychain",
        "/System/Library/Keychains/SystemRootCertificates.keychain",
    ]
    pem_parts = []
    for kc in keychains:
        if os.path.exists(kc):
            r = subprocess.run(
                ["security", "find-certificate", "-a", "-p", kc],
                capture_output=True,
                text=True,
            )
            if r.returncode == 0:
                pem_parts.append(r.stdout)

    # Append certifi's own bundle
    with open(certifi.where()) as f:
        pem_parts.append(f.read())

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".pem", delete=False, prefix="rcsb_ca_")
    tmp.write("\n".join(pem_parts))
    tmp.close()
    return tmp.name


_CA_BUNDLE = _build_ca_bundle()

_PDB_RE = __import__("re").compile(r"^[A-Z0-9]{4}$")


def _normalise_pdb_id(raw: str) -> str:
    """Recover PDB IDs mangled by Excel's auto-formatting.

    Excel reformats e.g. '6BHD' as '6,000 BHD' (thousands separator + space).
    Strip commas and collapse internal whitespace, then take the first
    four non-space characters that match the PDB pattern.
    """
    cleaned = raw.replace(",", "").replace(" ", "").strip()
    return cleaned[:4] if len(cleaned) >= 4 else cleaned


def _is_valid_pdb_id(pdb_id: str) -> bool:
    return bool(_PDB_RE.match(pdb_id.upper()))


# ---------------------------------------------------------------------------
# HTTP client
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Quality metrics
# ---------------------------------------------------------------------------


def get_entry_quality(client: RCSBClient, pdb_id: str) -> dict:
    data = client.get(f"{RCSB_DATA}/entry/{pdb_id}")
    if not data:
        return {}

    out = {}

    # Experimental method
    exptl = data.get("exptl") or []
    out["exp_method"] = exptl[0].get("method") if exptl else None

    # Resolution
    entry_info = data.get("rcsb_entry_info") or {}
    dres = entry_info.get("diffrn_resolution_high") or {}
    out["resolution_A"] = dres.get("value")

    # Refinement statistics
    refine = data.get("refine") or []
    if refine:
        out["r_work"] = refine[0].get("ls_R_factor_R_work")
        out["r_free"] = refine[0].get("ls_R_factor_R_free")
    else:
        out["r_work"] = None
        out["r_free"] = None

    # Geometry validation (rotamer/ramachandran/clashscore)
    geom = data.get("pdbx_vrpt_summary_geometry") or []
    if geom:
        g = geom[0]
        out["clashscore"] = g.get("clashscore")
        out["ramachandran_outliers_pct"] = g.get("percent_ramachandran_outliers")
        out["rotamer_outliers_pct"] = g.get("percent_rotamer_outliers")
        out["bonds_rmsz"] = g.get("bonds_RMSZ")
        out["angles_rmsz"] = g.get("angles_RMSZ")
    else:
        out["clashscore"] = None
        out["ramachandran_outliers_pct"] = None
        out["rotamer_outliers_pct"] = None
        out["bonds_rmsz"] = None
        out["angles_rmsz"] = None

    # Diffraction validation (RSRZ — X-ray only, requires EDS processing)
    diff = data.get("pdbx_vrpt_summary_diffraction") or []
    out["rsrz_outliers_pct"] = diff[0].get("percent_RSRZ_outliers") if diff else None

    # Non-polymer ligands
    ligands = entry_info.get("nonpolymer_bound_components") or []
    out["ligands_present"] = ",".join(ligands) if ligands else ""

    # Entity IDs (used later)
    identifiers = data.get("rcsb_entry_container_identifiers") or {}
    out["_polymer_entity_ids"] = identifiers.get("polymer_entity_ids") or []
    out["_nonpolymer_entity_ids"] = identifiers.get("non_polymer_entity_ids") or []

    return out


# ---------------------------------------------------------------------------
# Ligand interest classification
# ---------------------------------------------------------------------------

# CCD codes that are never drug-like ligands of interest.
# Ions, solvents, cryoprotectants, detergents, common cofactors that are
# not themselves the target of drug discovery.
_NON_INTERESTING_CCD: frozenset = frozenset({
    # Water / solvents / cryoprotectants
    "HOH", "DOD", "EDO", "EGL", "PEG", "PE3", "PE4", "PE5", "PE7", "PE8",
    "GOL", "DMS", "ACT", "ACE", "MPD", "IPA", "EOH", "PGE", "PGO", "TFP",
    "TRS", "BU3", "BU2", "BU1", "1BO", "BME", "2ME", "MLI",
    # Inorganic ions
    "MG", "ZN", "CA", "NA", "CL", "FE", "MN", "CO", "NI", "CU", "K",
    "BR", "IOD", "F", "CD", "HG", "PT", "AU", "AG", "LI", "RB", "CS",
    "BA", "SR", "AL", "GA", "IN", "TL", "PB", "BI", "SB", "AS", "SE",
    "TE", "PO4", "SO4", "SO3", "NO3", "NO2", "CO3", "CO2", "OH", "OXY",
    "FE2", "FE3", "FES", "SF4", "F3S", "CLF", "ZNO", "CUA", "CUB", "MO",
    "MO3", "MO5", "MO6", "W", "WO4", "VO4", "REO",
    # Detergents / lipids / amphiphiles
    "OLC", "OLA", "LMT", "BOG", "DDM", "DM", "OG", "NG", "HG", "C8E",
    "PLC", "LPC", "LPG", "LPE", "LPS", "PGP", "DAG", "TAG", "MAG",
    "SDS", "DPC", "LDAO", "CHAPSO", "CHAPS",
    # Buffers / additives
    "MES", "HEPES", "TRIS", "PIPES", "MOPS", "BICINE", "CAPS", "CHES",
    "EPPS", "BIS", "TCEP", "DTT", "BME", "EDO", "PEG", "MPD", "PG4",
    # Common non-drug cofactors / metals in enzyme active sites
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "CTP", "CDP", "UTP", "UDP",
    "NAD", "NAP", "NDP", "FAD", "FMN", "HEM", "HEC", "HEA", "COA", "COB",
    "SAM", "SAH", "THF", "TPP", "PLP", "PLR", "BIO", "FES", "SF4",
    "MTE", "F43", "CLA", "BCL", "PLM", "RET", "SRT",
})


def _is_interesting_ligand(comp_id: str, is_subject_of_investigation) -> bool:
    """Return True if the ligand is likely a drug-like molecule of interest.

    Uses the RCSB ISI (is_subject_of_investigation) flag as primary signal;
    falls back to CCD exclusion list for older entries where the flag is absent.
    """
    if is_subject_of_investigation is True:
        return True
    if is_subject_of_investigation is False:
        return False
    # Flag absent — use exclusion list
    return comp_id.upper() not in _NON_INTERESTING_CCD


# ---------------------------------------------------------------------------
# Ligand quality metrics + binding-site residue validation
# ---------------------------------------------------------------------------

# Sparse feature types: one element per outlier residue, values[0] = count
_SPARSE_FEATURE_TYPES = {
    "ANGLE_OUTLIERS", "BOND_OUTLIERS", "CLASHES", "SYMM_CLASHES",
    "RAMACHANDRAN_OUTLIER", "ROTAMER_OUTLIER", "CHIRAL_OUTLIERS",
    "PLANE_OUTLIERS", "STEREO_OUTLIER", "CIS-PEPTIDE",
    "MOGUL_ANGLE_OUTLIERS", "MOGUL_BOND_OUTLIERS",
}

# Feature types where values[] is a dense per-residue array starting at beg_seq_id
_DENSE_FEATURE_TYPES = {"RSRZ", "RSR", "RSCC", "OWAB", "Q_SCORE", "ASA"}


def _parse_residue_features(features: list) -> dict:
    """Return {seq_id: {feature_type: value}} from rcsb_polymer_instance_feature."""
    residue_data: dict = {}
    for feat in features:
        ftype = feat.get("type", "")
        for pos in feat.get("feature_positions") or []:
            values = pos.get("values") or []
            beg = pos.get("beg_seq_id")
            if beg is None or not values:
                continue
            if ftype in _DENSE_FEATURE_TYPES:
                for i, val in enumerate(values):
                    rid = beg + i
                    residue_data.setdefault(rid, {})[ftype] = val
            elif ftype in _SPARSE_FEATURE_TYPES:
                residue_data.setdefault(beg, {})[ftype] = values[0]
    return residue_data


def _traffic_light(
    rscc,
    rsr,
    rmsz_bonds,
    rmsz_angles,
    intermolecular_clashes,
    contact_outlier_fraction: float,
) -> str:
    """Aggregate per-ligand and per-contact-residue metrics into a traffic-light score.

    Scoring table (each criterion contributes 0=good, 1=fair, 2=bad):
      RSCC       ≥0.80→0   ≥0.60→1   <0.60→2   (hard override: <0.50 → "bad")
      RSR        ≤0.20→0   ≤0.35→1   >0.35→2
      RMSZ bonds ≤1.5→0    ≤2.5→1    >2.5→2
      RMSZ angles≤1.5→0    ≤2.5→1    >2.5→2
      i-mol clashes 0→0    1-2→1     ≥3→2
      contact outlier fraction ≤0.10→0 ≤0.25→1 >0.25→2

    Mean score: <0.67→"good"  <1.33→"fair"  ≥1.33→"bad"
    Any missing metric is simply omitted from the average.
    """
    # Hard override
    if rscc is not None and rscc < 0.50:
        return "bad"

    scores = []
    if rscc is not None:
        scores.append(0 if rscc >= 0.80 else (1 if rscc >= 0.60 else 2))
    if rsr is not None:
        scores.append(0 if rsr <= 0.20 else (1 if rsr <= 0.35 else 2))
    if rmsz_bonds is not None:
        scores.append(0 if rmsz_bonds <= 1.5 else (1 if rmsz_bonds <= 2.5 else 2))
    if rmsz_angles is not None:
        scores.append(0 if rmsz_angles <= 1.5 else (1 if rmsz_angles <= 2.5 else 2))
    if intermolecular_clashes is not None:
        scores.append(0 if intermolecular_clashes == 0 else (1 if intermolecular_clashes <= 2 else 2))
    scores.append(0 if contact_outlier_fraction <= 0.10 else (1 if contact_outlier_fraction <= 0.25 else 2))

    if not scores:
        return ""
    mean = sum(scores) / len(scores)
    return "good" if mean < 0.67 else ("fair" if mean < 1.33 else "bad")


def _get_contact_residue_outlier_fraction(
    client: RCSBClient,
    pdb_id: str,
    neighbors: list,
    cache: dict,
) -> float:
    """Fraction of contact residues (≤CONTACT_DIST_CUTOFF Å) that carry any validation outlier.

    `cache` is a per-entry dict keyed by asym_id to avoid re-fetching the same chain.
    An outlier is any entry in _SPARSE_FEATURE_TYPES with a non-zero value.
    """
    contact_seq_ids: dict = {}  # asym_id → set of seq_ids within cutoff
    for n in neighbors:
        if (n.get("distance") or 999) > CONTACT_DIST_CUTOFF:
            continue
        asym = n.get("target_asym_id")
        seq_id = n.get("target_seq_id")
        if asym and seq_id is not None:
            contact_seq_ids.setdefault(asym, set()).add(seq_id)

    if not contact_seq_ids:
        return 0.0

    total = 0
    outliers = 0
    for asym, seq_ids in contact_seq_ids.items():
        if asym not in cache:
            inst = client.get(f"{RCSB_DATA}/polymer_entity_instance/{pdb_id}/{asym}")
            features = (inst.get("rcsb_polymer_instance_feature") or []) if inst else []
            cache[asym] = _parse_residue_features(features)
        res_data = cache[asym]
        for sid in seq_ids:
            total += 1
            rdata = res_data.get(sid, {})
            if any(k in _SPARSE_FEATURE_TYPES for k in rdata):
                outliers += 1

    return outliers / total if total else 0.0


def get_ligand_quality(client: RCSBClient, pdb_id: str, nonpolymer_entity_ids: list) -> list:
    """Return per-ligand quality metrics including binding-site residue validation.

    For each ligand instance:
    - Fetches rcsb_nonpolymer_instance_validation_score (RSCC, RSR, Mogul RMSZ, clashes)
    - Fetches rcsb_target_neighbors and identifies contact residues ≤ CONTACT_DIST_CUTOFF Å
    - Fetches polymer_entity_instance feature data for each contact chain (cached)
    - Computes the fraction of contact residues carrying any geometry/density outlier
    - Derives a single traffic-light score: "good", "fair", or "bad"
    """
    results = []
    chain_feature_cache: dict = {}  # asym_id → parsed residue feature dict, shared across ligands

    for eid in nonpolymer_entity_ids:
        entity_data = client.get(f"{RCSB_DATA}/nonpolymer_entity/{pdb_id}/{eid}")
        if not entity_data:
            continue
        container = entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}
        comp_id = container.get("comp_id") or eid
        instance_ids = container.get("nonpolymer_entity_instance_ids") or []

        for asym_id in instance_ids:
            inst = client.get(f"{RCSB_DATA}/nonpolymer_entity_instance/{pdb_id}/{asym_id}")
            if not inst:
                continue

            scores = (inst.get("rcsb_nonpolymer_instance_validation_score") or [None])[0] or {}
            rscc = scores.get("RSCC")
            rsr = scores.get("RSR")
            rmsz_bonds = scores.get("mogul_bonds_RMSZ")
            rmsz_angles = scores.get("mogul_angles_RMSZ")
            intermolecular_clashes = scores.get("intermolecular_clashes")
            isi = scores.get("is_subject_of_investigation")  # True/False/None
            interesting = _is_interesting_ligand(comp_id, isi)

            neighbors = inst.get("rcsb_target_neighbors") or []
            outlier_frac = _get_contact_residue_outlier_fraction(
                client, pdb_id, neighbors, chain_feature_cache
            )

            n_contacts = sum(
                1 for n in neighbors
                if (n.get("distance") or 999) <= CONTACT_DIST_CUTOFF
            )
            contact_residues = sorted({
                f"{n.get('target_comp_id','?')}{n.get('target_auth_seq_id','?')}"
                f"({n.get('target_asym_id','')})"
                for n in neighbors
                if (n.get("distance") or 999) <= CONTACT_DIST_CUTOFF
            })

            tl = _traffic_light(
                rscc, rsr, rmsz_bonds, rmsz_angles,
                intermolecular_clashes, outlier_frac,
            )

            results.append({
                "ligand_id": comp_id,
                "chain_id": asym_id,
                "is_interesting": interesting,
                "ligand_rscc": rscc,
                "ligand_rsr": rsr,
                "ligand_rmsz_bonds": rmsz_bonds,
                "ligand_rmsz_angles": rmsz_angles,
                "ligand_intermolecular_clashes": intermolecular_clashes,
                "contact_residue_count": n_contacts,
                "contact_outlier_fraction": round(outlier_frac, 3),
                "contact_residues": ";".join(contact_residues),
                "binding_quality": tl,
            })

    return results


# ---------------------------------------------------------------------------
# Polymer entity data: UniProt IDs, sequence, cofactors, binding features
# ---------------------------------------------------------------------------


def get_polymer_entities(client: RCSBClient, pdb_id: str, entity_ids: list[str]) -> list[dict]:
    entities = []
    for eid in entity_ids:
        data = client.get(f"{RCSB_DATA}/polymer_entity/{pdb_id}/{eid}")
        if not data:
            continue

        entity: dict = {}

        # Polymer type filter — only protein chains
        entity_poly = data.get("entity_poly") or {}
        if entity_poly.get("rcsb_entity_polymer_type") != "Protein":
            continue

        # UniProt IDs
        container = data.get("rcsb_polymer_entity_container_identifiers") or {}
        entity["uniprot_ids"] = container.get("uniprot_ids") or []

        # Canonical sequence
        entity["sequence"] = entity_poly.get("pdbx_seq_one_letter_code_can") or ""

        # Direct binders from ChEMBL/DrugBank (neighbor_flag=="N" means this protein is the target)
        cofactors = []
        for c in data.get("rcsb_target_cofactors") or []:
            if c.get("neighbor_flag") == "Y":
                continue
            cofactors.append(
                {
                    "name": c.get("cofactor_name"),
                    "smiles": c.get("cofactor_SMILES"),
                    "inchikey": c.get("cofactor_InChIKey"),
                    "chem_comp_id": c.get("cofactor_chem_comp_id"),
                    "source": c.get("resource_name", ""),
                }
            )
        entity["cofactors"] = cofactors

        # Binding site / functional site features from UniProt/CSA
        site_types = {"BINDING_SITE", "ACTIVE_SITE", "METAL_COORDINATION", "SITE"}
        site_features = []
        for f in data.get("rcsb_polymer_entity_feature") or []:
            if f.get("type") in site_types:
                site_features.append(
                    {
                        "name": f.get("name"),
                        "feature_id": f.get("feature_id"),
                        "source": f.get("provenance_source"),
                        "type": f.get("type"),
                    }
                )
        entity["site_features"] = site_features

        entities.append(entity)

    return entities


# ---------------------------------------------------------------------------
# Binder lookup for related entries
# ---------------------------------------------------------------------------

_MAX_RELATED_BINDER_ENTRIES = 5  # cap API calls for sibling/full-length binder lookups


def _extract_direct_binders(client: RCSBClient, pdb_id: str) -> list[dict]:
    """Return all direct (neighbor_flag=='N') binders from all protein entities of pdb_id."""
    entry = client.get(f"{RCSB_DATA}/entry/{pdb_id}")
    if not entry:
        return []
    entity_ids = (
        (entry.get("rcsb_entry_container_identifiers") or {}).get("polymer_entity_ids") or []
    )
    binders = []
    for eid in entity_ids:
        data = client.get(f"{RCSB_DATA}/polymer_entity/{pdb_id}/{eid}")
        if not data:
            continue
        if (data.get("entity_poly") or {}).get("rcsb_entity_polymer_type") != "Protein":
            continue
        for c in data.get("rcsb_target_cofactors") or []:
            if c.get("neighbor_flag") == "Y":
                continue
            binders.append(
                {
                    "name": c.get("cofactor_name"),
                    "smiles": c.get("cofactor_SMILES"),
                    "inchikey": c.get("cofactor_InChIKey"),
                    "chem_comp_id": c.get("cofactor_chem_comp_id"),
                    "source": c.get("resource_name", ""),
                    "from_pdb": pdb_id,
                }
            )
    return binders


# ---------------------------------------------------------------------------
# PDBe binding sites (SIFTS-curated)
# ---------------------------------------------------------------------------


def get_pdbe_binding_sites(client: RCSBClient, pdb_id: str) -> list[dict]:
    data = client.get(f"{PDBE_API}/{pdb_id.lower()}")
    if not data:
        return []

    entry_data = data.get(pdb_id.lower()) or []
    sites = []
    for site in entry_data:
        sites.append(
            {
                "description": site.get("site_description") or site.get("details", ""),
                "source": "PDBe-KB/SIFTS",
            }
        )
    return sites


# ---------------------------------------------------------------------------
# Holo structure lookup for known binders
# ---------------------------------------------------------------------------

_MAX_HOLO_ENTRIES = 3   # holo structures evaluated per binder (API cost control)


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
) -> list[str]:
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
    binders: list[dict],
    uniprot_ids: list[str],
) -> list[dict]:
    """For each known binder, find holo PDB entries and score ligand binding quality.

    Steps per binder:
    1. Resolve InChIKey → CCD code (skip if no InChIKey and no chem_comp_id)
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
            # Keep only instances whose CCD matches this binder
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


# ---------------------------------------------------------------------------
# Related PDB entries
# ---------------------------------------------------------------------------

# Chains longer than query × this factor are treated as "full-length" proteins.
_SEQ_LEN_RATIO = 1.4


def _uniprot_length_search(
    client: RCSBClient,
    uniprot_id: str,
    length_operator: str,
    length_value: int,
    max_rows: int,
) -> list[str]:
    """UniProt exact-match AND entity sequence-length comparison."""
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_id,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": length_operator,
                        "value": length_value,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows}},
    }
    data = client.post(RCSB_SEARCH, payload)
    if not data:
        return []
    return [r["identifier"] for r in data.get("result_set") or []]


def get_related_by_uniprot(client: RCSBClient, uniprot_id: str, max_rows: int = 200) -> list[str]:
    """All entries sharing uniprot_id (no length split — used when query seq len is unknown)."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows}},
    }
    data = client.post(RCSB_SEARCH, payload)
    if not data:
        return []
    return [r["identifier"] for r in data.get("result_set") or []]


def get_related_by_uniprot_split(
    client: RCSBClient,
    uniprot_id: str,
    query_seq_len: int,
    max_rows: int = 200,
) -> tuple[list[str], list[str]]:
    """Split same-UniProt entries into (siblings, full_length) by chain length.

    'Siblings' share the same region (chain length ≤ query × _SEQ_LEN_RATIO).
    'Full-length' have significantly longer chains (> that threshold).
    Two search requests; no per-hit extra fetches.
    """
    threshold = int(query_seq_len * _SEQ_LEN_RATIO)
    siblings = _uniprot_length_search(client, uniprot_id, "less_or_equal", threshold, max_rows)
    full_length = _uniprot_length_search(client, uniprot_id, "greater", threshold, max_rows)
    return siblings, full_length


def get_related_by_sequence(
    client: RCSBClient,
    sequence: str,
    identity_cutoff: float = 0.9,
    max_rows: int = 25,
) -> list[str]:
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.1,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": sequence,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows}},
    }
    data = client.post(RCSB_SEARCH, payload)
    if not data:
        return []
    return [r["identifier"] for r in data.get("result_set") or []]


# ---------------------------------------------------------------------------
# Row enrichment
# ---------------------------------------------------------------------------


def enrich_row(
    row: dict,
    client: RCSBClient,
    pdb_col: str,
    uniprot_col: str | None,
    seq_identity: float,
    max_related: int,
) -> dict:
    pdb_id = _normalise_pdb_id(str(row[pdb_col]).strip()).upper()
    result = dict(row)

    valid_pdb = _is_valid_pdb_id(pdb_id)
    if not valid_pdb:
        log.warning("[%s] Not a standard PDB ID — skipping structure lookups", pdb_id)

    # --- Quality metrics (only for valid PDB IDs) ---
    quality: dict = {}
    entities: list = []

    ligand_metrics: list = []

    if valid_pdb:
        log.info("[%s] Fetching quality metrics", pdb_id)
        quality = get_entry_quality(client, pdb_id)
        entity_ids = quality.pop("_polymer_entity_ids", [])
        nonpolymer_entity_ids = quality.pop("_nonpolymer_entity_ids", [])

        log.info("[%s] Fetching %d polymer entities", pdb_id, len(entity_ids))
        entities = get_polymer_entities(client, pdb_id, entity_ids)

        if nonpolymer_entity_ids:
            log.info("[%s] Fetching ligand quality for %d entities", pdb_id, len(nonpolymer_entity_ids))
            ligand_metrics = get_ligand_quality(client, pdb_id, nonpolymer_entity_ids)

    for col in (
        "exp_method",
        "resolution_A",
        "r_work",
        "r_free",
        "clashscore",
        "ramachandran_outliers_pct",
        "rotamer_outliers_pct",
        "rsrz_outliers_pct",
        "bonds_rmsz",
        "angles_rmsz",
        "ligands_present",
    ):
        result[col] = quality.get(col)

    # Split ligand metrics into interesting (drug-like) vs non-interesting (ions/solvents/cofactors)
    interesting_ligands = [m for m in ligand_metrics if m.get("is_interesting")]
    noninteresting_ligands = [m for m in ligand_metrics if not m.get("is_interesting")]

    result["ligand_quality"] = json.dumps(interesting_ligands) if interesting_ligands else ""
    result["cofactor_ion_quality"] = json.dumps(noninteresting_ligands) if noninteresting_ligands else ""

    # Convenience summary columns from interesting ligands only
    result["ligands_interesting"] = ",".join(
        dict.fromkeys(m["ligand_id"] for m in interesting_ligands)
    )
    result["ligands_noninteresting"] = ",".join(
        dict.fromkeys(m["ligand_id"] for m in noninteresting_ligands)
    )

    # Collect all UniProt IDs: from input column first, then from resolved entities
    all_uniprot: list[str] = []
    if uniprot_col and row.get(uniprot_col):
        all_uniprot = [str(row[uniprot_col]).strip()]
    for e in entities:
        for uid in e["uniprot_ids"]:
            if uid not in all_uniprot:
                all_uniprot.append(uid)

    # Back-fill UniProt: write into existing column if empty, or add resolved_uniprot column
    if all_uniprot:
        if uniprot_col and not row.get(uniprot_col):
            result[uniprot_col] = ",".join(all_uniprot)
            log.info("[%s] Resolved UniProt from PDB entry: %s", pdb_id, result[uniprot_col])
        elif not uniprot_col:
            result["resolved_uniprot"] = ",".join(all_uniprot)

    # Collect sequence from first entity (for fallback sequence search)
    first_sequence = next((e["sequence"] for e in entities if e.get("sequence")), "")

    # --- Related entries ---
    sibling_ids: list[str] = []
    fulllength_ids: list[str] = []
    search_method = ""

    if all_uniprot:
        log.info("[%s] Searching related entries by UniProt %s", pdb_id, all_uniprot[0])
        if first_sequence:
            sibling_ids, fulllength_ids = get_related_by_uniprot_split(
                client, all_uniprot[0], len(first_sequence), max_rows=max_related * 2
            )
            search_method = "uniprot_split"
        else:
            # No sequence available (e.g. UniProt from input column only) — no split possible
            sibling_ids = get_related_by_uniprot(client, all_uniprot[0], max_rows=max_related * 2)
            search_method = "uniprot"
    elif first_sequence:
        log.info("[%s] No UniProt; falling back to sequence similarity search", pdb_id)
        sibling_ids = get_related_by_sequence(
            client, first_sequence, identity_cutoff=seq_identity, max_rows=max_related
        )
        search_method = f"sequence_id_{seq_identity}"
    else:
        search_method = "none"

    # Exclude the input PDB ID itself from both lists
    sibling_ids = [r for r in sibling_ids if r.upper() != pdb_id]
    fulllength_ids = [r for r in fulllength_ids if r.upper() != pdb_id]

    # Rule: when only one list is non-empty, use related_pdb_ids (column 1)
    if sibling_ids or not fulllength_ids:
        result["related_pdb_ids"] = ",".join(sibling_ids[:max_related])
        result["related_pdb_count"] = len(sibling_ids)
    else:
        result["related_pdb_ids"] = ""
        result["related_pdb_count"] = 0

    result["fulllength_pdb_ids"] = ",".join(fulllength_ids[:max_related]) if fulllength_ids else ""
    result["fulllength_pdb_count"] = len(fulllength_ids)
    result["related_search_method"] = search_method

    # --- Binding sites ---
    all_site_features = []
    binding_sources: list[str] = []

    for e in entities:
        all_site_features.extend(e.get("site_features") or [])

    pdbe_sites = get_pdbe_binding_sites(client, pdb_id) if valid_pdb else []

    # --- Known binders (direct targets only, no neighbour annotations) ---
    # CCD codes co-crystallised in the query structure — used as fragment-coverage evidence
    query_ccd_codes = set(quality.get("ligands_present", "").split(",")) - {""}

    # 1. Direct binders from the query structure's own polymer entities
    all_cofactors: list[dict] = []
    for e in entities:
        all_cofactors.extend(e.get("cofactors") or [])

    # 2. Direct binders from sibling structures (same region — trust unconditionally)
    for sid in sibling_ids[:_MAX_RELATED_BINDER_ENTRIES]:
        log.info("[%s] Fetching binders from sibling %s", pdb_id, sid)
        for b in _extract_direct_binders(client, sid):
            b["binder_source_type"] = "sibling"
            all_cofactors.append(b)

    # 3. Direct binders from full-length structures — only accept if the binder's CCD code
    #    matches a ligand co-crystallised in the query fragment (evidence it binds that domain)
    for fid in fulllength_ids[:_MAX_RELATED_BINDER_ENTRIES]:
        log.info("[%s] Fetching binders from full-length %s", pdb_id, fid)
        for b in _extract_direct_binders(client, fid):
            if b.get("chem_comp_id") and b["chem_comp_id"] in query_ccd_codes:
                b["binder_source_type"] = "fulllength"
                all_cofactors.append(b)

    # Deduplicate by InChIKey (fall back to name when InChIKey absent)
    seen: set[str] = set()
    unique_cofactors: list[dict] = []
    for c in all_cofactors:
        key = c.get("inchikey") or c.get("name") or ""
        if key and key not in seen:
            seen.add(key)
            unique_cofactors.append(c)

    if all_site_features:
        binding_sources.append("UniProt/CSA")
    if unique_cofactors:
        binding_sources.append("ChEMBL/DrugBank")
    if pdbe_sites:
        binding_sources.append("PDBe-KB/SIFTS")

    result["has_binding_site"] = bool(binding_sources)
    result["binding_site_sources"] = ",".join(binding_sources)

    site_notes_parts = []
    if all_site_features:
        names = [s["name"] for s in all_site_features if s.get("name")]
        if names:
            site_notes_parts.append("; ".join(dict.fromkeys(names)))
    if pdbe_sites:
        descs = [s["description"] for s in pdbe_sites if s.get("description")]
        if descs:
            site_notes_parts.append("; ".join(descs[:3]))
    result["binding_site_notes"] = " | ".join(site_notes_parts) if site_notes_parts else ""

    if unique_cofactors:
        result["known_binders"] = ",".join(
            c["name"] for c in unique_cofactors if c.get("name")
        )
        result["known_binder_smiles"] = ",".join(
            c["smiles"] for c in unique_cofactors if c.get("smiles")
        )
    else:
        result["known_binders"] = ""
        result["known_binder_smiles"] = ""

    # --- Holo structure binding quality ---
    # For each known binder, find co-crystallised PDB entries and score quality.
    # This works even when the query structure is apo or a fragment with no ligand.
    _rank = {"bad": 2, "fair": 1, "good": 0}

    holo_results: list[dict] = []
    if unique_cofactors and all_uniprot:
        log.info("[%s] Looking up holo structures for %d known binder(s)", pdb_id, len(unique_cofactors))
        holo_results = get_holo_ligand_quality(client, unique_cofactors, all_uniprot)

    # Restrict traffic-light aggregation to interesting ligands only
    query_tl_values = [m["binding_quality"] for m in interesting_ligands if m.get("binding_quality")]
    holo_tl_values = [
        m["binding_quality"]
        for hr in holo_results
        for m in hr.get("holo_ligand_metrics", [])
        if m.get("binding_quality")
    ]
    all_tl = query_tl_values + holo_tl_values
    result["ligand_binding_quality"] = (
        max(all_tl, key=lambda x: _rank.get(x, -1)) if all_tl else ""
    )

    # Holo PDB entries keyed by binder name, serialised as JSON
    result["holo_quality"] = json.dumps(holo_results) if holo_results else ""

    return result


# ---------------------------------------------------------------------------
# CSV column detection
# ---------------------------------------------------------------------------


def detect_pdb_col(columns: list[str], hint: str | None) -> str:
    if hint:
        if hint in columns:
            return hint
        raise ValueError(f"Column '{hint}' not found. Available: {columns}")
    for col in columns:
        if col.strip().lower() in _PDB_COL_ALIASES:
            return col
    raise ValueError(f"Cannot auto-detect PDB ID column. Use --pdb-col. Available: {columns}")


def detect_uniprot_col(columns: list[str], hint: str | None) -> str | None:
    if hint:
        if hint in columns:
            return hint
        log.warning("UniProt column '%s' not found; skipping UniProt lookup", hint)
        return None
    for col in columns:
        if col.strip().lower() in {"uniprot", "uniprot_id", "uniprot_acc", "accession"}:
            return col
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Enrich a CSV of PDB IDs with RCSB quality metrics and binding site annotations."
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV file path")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")
    parser.add_argument("-p", "--pdb-col", default=None, help="Column name containing PDB IDs")
    parser.add_argument(
        "-u",
        "--uniprot-col",
        default=None,
        help="Column name containing UniProt accessions (optional)",
    )
    parser.add_argument(
        "--seq-identity",
        type=float,
        default=0.9,
        help="Sequence identity cutoff for fallback search (default: 0.9)",
    )
    parser.add_argument(
        "--max-related",
        type=int,
        default=25,
        help="Max related PDB entries to return (default: 25)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.1,
        help="Seconds to wait between API requests (default: 0.1)",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    pdb_col = detect_pdb_col(df.columns.tolist(), args.pdb_col)
    uniprot_col = detect_uniprot_col(df.columns.tolist(), args.uniprot_col)

    log.info(
        "Input: %d rows | PDB column: '%s' | UniProt column: %s",
        len(df),
        pdb_col,
        uniprot_col or "none",
    )

    client = RCSBClient(delay=args.delay)
    enriched_rows = []

    for i, row in enumerate(df.to_dict(orient="records"), 1):
        pdb_id = _normalise_pdb_id(str(row.get(pdb_col, "")).strip()).upper()
        if not pdb_id:
            log.warning("Row %d: empty PDB ID, skipping", i)
            enriched_rows.append(row)
            continue
        if pdb_id != str(row.get(pdb_col, "")).strip().upper():
            log.warning("Row %d: normalised PDB ID '%s' → '%s' (Excel formatting?)",
                        i, row.get(pdb_col, ""), pdb_id)
        log.info("Processing row %d/%d: %s", i, len(df), pdb_id)
        try:
            enriched = enrich_row(
                row=row,
                client=client,
                pdb_col=pdb_col,
                uniprot_col=uniprot_col,
                seq_identity=args.seq_identity,
                max_related=args.max_related,
            )
        except Exception as exc:
            log.error("[%s] Failed: %s", pdb_id, exc, exc_info=True)
            enriched = dict(row)
            for col in (
                "exp_method",
                "resolution_A",
                "r_work",
                "r_free",
                "clashscore",
                "ramachandran_outliers_pct",
                "rotamer_outliers_pct",
                "rsrz_outliers_pct",
                "bonds_rmsz",
                "angles_rmsz",
                "ligands_present",
                "ligands_interesting",
                "ligands_noninteresting",
                "ligand_quality",
                "cofactor_ion_quality",
                "ligand_binding_quality",
                "holo_quality",
                "related_pdb_ids",
                "related_pdb_count",
                "fulllength_pdb_ids",
                "fulllength_pdb_count",
                "related_search_method",
                "has_binding_site",
                "binding_site_sources",
                "binding_site_notes",
                "known_binders",
                "known_binder_smiles",
            ):
                enriched.setdefault(col, None)
        enriched_rows.append(enriched)

    out_df = pd.DataFrame(enriched_rows)
    out_df.to_csv(args.output, index=False)
    log.info("Written %d rows to %s", len(out_df), args.output)


if __name__ == "__main__":
    main()
