"""Row-level orchestration: fetch all data for a single PDB ID and assemble output columns."""

import json
import logging
import re

from .binding_sites import get_pdbe_binding_sites
from .client import RCSBClient
from .entities import extract_direct_binders, get_polymer_entities
from .holo import get_holo_ligand_quality
from .quality import get_entry_quality, get_ligand_quality, iridium_score
from .related import (
    get_related_by_sequence,
    get_related_by_uniprot,
    get_related_by_uniprot_split,
)

log = logging.getLogger(__name__)

_MAX_RELATED_BINDER_ENTRIES = 5  # cap API calls for sibling/full-length binder lookups
_SEARCH_MAX_ROWS = 1000           # upper bound passed to RCSB search; independent of max_related


def _entity_matches_names(description: str, name_filters: list) -> bool:
    """Return True if any filter term appears as a whole whitespace-delimited token
    in *description* (case-insensitive).  Hyphenated words like 'Tubulin-Tyrosine'
    are a single token and will NOT match 'Tubulin'.
    """
    tokens = set(description.lower().split())
    return any(term.lower() in tokens for term in name_filters)


def _collect_entity_names(entities: list, peptide_entities: list, ligand_metrics: list) -> str:
    """Build a comma-separated list of all entity descriptions for a PDB entry.

    Polymer descriptions come from entity['description'] (set by get_polymer_entities).
    Non-polymer descriptions come from ligand_metric['description'] (set by get_ligand_quality).
    Only unique, non-empty names are included; order is polymer entities first,
    then non-polymer entities.
    """
    seen: set = set()
    names: list = []
    for e in list(entities) + list(peptide_entities):
        d = e.get("description") or ""
        if d and d not in seen:
            seen.add(d)
            names.append(d)
    for m in ligand_metrics:
        d = m.get("description") or ""
        if d and d not in seen:
            seen.add(d)
            names.append(d)
    return ",".join(names)


_PDB_RE = re.compile(r"^[A-Z0-9]{4}$")

# Columns that belong to ligand sub-rows only.  Primary rows leave these None.
LIGAND_DETAIL_COLS = (
    "ligand_type",
    "ligand_id",
    "chain_id",
    "ligand_rscc",
    "ligand_rsr",
    "ligand_rmsz_bonds",
    "ligand_rmsz_angles",
    "ligand_intermolecular_clashes",
    "contact_residue_count",
    "contact_outlier_fraction",
    "contact_residues",
    "binding_quality",
)


def _normalise_pdb_id(raw: str) -> str:
    """Recover PDB IDs mangled by Excel's auto-formatting.

    Excel reformats e.g. '6BHD' as '6,000 BHD' (thousands separator + space).
    Strip commas and collapse internal whitespace, then take the first four characters.
    """
    cleaned = raw.replace(",", "").replace(" ", "").strip()
    return cleaned[:4] if len(cleaned) >= 4 else cleaned


def _is_valid_pdb_id(pdb_id: str) -> bool:
    return bool(_PDB_RE.match(pdb_id.upper()))


_CRYSTAL_QUALITY_COLS = (
    "exp_method", "resolution_A", "r_work", "r_free", "clashscore",
    "ramachandran_outliers_pct", "rotamer_outliers_pct", "rsrz_outliers_pct",
    "bonds_rmsz", "angles_rmsz", "ligands_present",
)


def _fetch_related_ligand_data(client: RCSBClient, pdb_id: str) -> dict:
    """Fetch ligand and crystal quality data for a single related PDB entry.

    Returns a dict with keys:
        pdb_id, entry_quality, ligand_metrics, peptide_entities, has_ligands
    A 404 response (quality={}) immediately returns has_ligands=False.
    """
    quality = get_entry_quality(client, pdb_id)
    if not quality:
        return {
            "pdb_id": pdb_id, "entry_quality": {},
            "ligand_metrics": [], "peptide_entities": [], "has_ligands": False,
        }

    entity_ids = quality.pop("_polymer_entity_ids", [])
    nonpolymer_entity_ids = quality.pop("_nonpolymer_entity_ids", [])
    entry_quality = {col: quality.get(col) for col in _CRYSTAL_QUALITY_COLS}

    ligand_metrics = []
    if nonpolymer_entity_ids:
        ligand_metrics = get_ligand_quality(client, pdb_id, nonpolymer_entity_ids)

    peptide_entities = []
    receptor_entities = []
    species = ""
    if entity_ids:
        all_entities = get_polymer_entities(client, pdb_id, entity_ids)
        peptide_entities = [e for e in all_entities if e.get("is_peptide_ligand")]
        receptor_entities = [e for e in all_entities if not e.get("is_peptide_ligand")]
        species_labels = [e["species"] for e in receptor_entities if e.get("species")]
        species = ",".join(dict.fromkeys(species_labels))

    # A ligand is meaningful by the same criterion used in build_ligand_rows
    interesting = [
        m for m in ligand_metrics
        if m.get("is_interesting") and (m.get("contact_residue_count") or m.get("binding_quality"))
    ]
    has_ligands = bool(interesting or peptide_entities)

    _tl_order = {"good": 0, "fair": 1, "bad": 2}
    _best_ligand_traffic = min(
        (m["binding_quality"] for m in interesting if m.get("binding_quality")),
        key=lambda v: _tl_order.get(v, 3),
        default=None,
    )

    entity_names = _collect_entity_names(receptor_entities, peptide_entities, ligand_metrics)
    _grade, _lig_used = iridium_score(entry_quality, _best_ligand_traffic)

    return {
        "pdb_id": pdb_id,
        "entry_quality": entry_quality,
        "ligand_metrics": ligand_metrics,
        "peptide_entities": peptide_entities,
        "has_ligands": has_ligands,
        "species": species,
        "entity_names": entity_names,
        "structure_quality": _grade,
        "structure_quality_ligand_used": _lig_used,
    }


def build_ligand_rows(
    pdb_id: str,
    ligand_metrics: list,
    peptide_entities: list,
    all_output_cols: tuple,
    tags: dict = None,
) -> list:
    """Return one sub-row dict per meaningfully-bound ligand/peptide.

    Small-molecule criterion: is_interesting AND (contact_residue_count > 0 OR binding_quality != "").
    Peptide ligands: always included (annotation itself implies binding).
    All columns not specific to ligands are set to None for visual clarity in the output CSV.
    `tags` is an optional dict of extra column values to set on every returned row (e.g. the
    related_pdb_ids or fulllength_pdb_ids tag).
    """
    rows = []
    blank = {col: None for col in all_output_cols}

    for m in ligand_metrics:
        if not m.get("is_interesting"):
            continue
        if not m.get("contact_residue_count") and not m.get("binding_quality"):
            continue
        row = dict(blank)
        row["row_type"] = "ligand"
        row["parent_pdb_id"] = pdb_id
        row["ligand_type"] = "small_molecule"
        for col in LIGAND_DETAIL_COLS:
            if col != "ligand_type":
                row[col] = m.get(col)
        if tags:
            row.update(tags)
        rows.append(row)

    for pe in peptide_entities:
        label = pe.get("bird_id") or pe.get("sequence") or ""
        if not label:
            continue
        row = dict(blank)
        row["row_type"] = "ligand"
        row["parent_pdb_id"] = pdb_id
        row["ligand_type"] = "peptide"
        row["ligand_id"] = label
        if tags:
            row.update(tags)
        rows.append(row)

    return rows


def build_related_row(
    pdb_id: str,
    entry_data: dict,
    all_output_cols: tuple,
    tags: dict = None,
) -> dict:
    """Return a single row for a related entry that has no qualifying ligand.

    Carries structure quality metrics and provenance tags only; all ligand-
    specific columns are None.  row_type is 'related'.
    """
    row = {col: None for col in all_output_cols}
    row["row_type"] = "related"
    row["parent_pdb_id"] = pdb_id
    row["structure_quality"] = entry_data.get("structure_quality")
    row["structure_quality_ligand_used"] = entry_data.get("structure_quality_ligand_used")
    row["species"] = entry_data.get("species")
    row["entity_names"] = entry_data.get("entity_names", "")
    row.update(entry_data.get("entry_quality", {}))
    if tags:
        row.update(tags)
    return row


def enrich_row(
    row: dict,
    client: RCSBClient,
    pdb_col: str,
    uniprot_col: str | None,
    seq_identity: float,
    max_related: int,
    input_pdb_ids: frozenset = frozenset(),
    entity_name_filters: list | None = None,
    related_data_cache: dict | None = None,
    binders_cache: dict | None = None,
    uniprot_search_cache: dict | None = None,
    include_all_related: bool = False,
) -> dict:
    pdb_id = _normalise_pdb_id(str(row[pdb_col]).strip()).upper()
    result = dict(row)

    valid_pdb = _is_valid_pdb_id(pdb_id)
    if not valid_pdb:
        log.warning("[%s] Not a standard PDB ID — skipping structure lookups", pdb_id)

    quality: dict = {}
    entities: list = []
    ligand_metrics: list = []
    peptide_entities: list = []

    if valid_pdb:
        log.info("[%s] Fetching quality metrics", pdb_id)
        quality = get_entry_quality(client, pdb_id)
        entity_ids = quality.pop("_polymer_entity_ids", [])
        nonpolymer_entity_ids = quality.pop("_nonpolymer_entity_ids", [])

        log.info("[%s] Fetching %d polymer entities", pdb_id, len(entity_ids))
        all_entities = get_polymer_entities(client, pdb_id, entity_ids)
        peptide_entities = [e for e in all_entities if e.get("is_peptide_ligand")]
        entities = [e for e in all_entities if not e.get("is_peptide_ligand")]

        if entity_name_filters:
            before = len(entities)
            entities = [
                e for e in entities
                if _entity_matches_names(e.get("description", ""), entity_name_filters)
            ]
            log.info(
                "[%s] Entity name filter %s: kept %d of %d receptor entities",
                pdb_id, entity_name_filters, len(entities), before,
            )

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

    interesting_ligands = [m for m in ligand_metrics if m.get("is_interesting")]
    noninteresting_ligands = [m for m in ligand_metrics if not m.get("is_interesting")]

    result["ligand_quality"] = json.dumps(interesting_ligands) if interesting_ligands else ""

    result["ligands_interesting"] = ",".join(
        dict.fromkeys(m["ligand_id"] for m in interesting_ligands)
    )
    result["ligands_noninteresting"] = ",".join(
        dict.fromkeys(m["ligand_id"] for m in noninteresting_ligands)
    )

    species_labels = [e["species"] for e in entities if e.get("species")]
    result["species"] = ",".join(dict.fromkeys(species_labels))

    result["entity_names"] = _collect_entity_names(entities, peptide_entities, ligand_metrics)

    _tl_order = {"good": 0, "fair": 1, "bad": 2}
    _best_ligand_traffic = min(
        (m["binding_quality"] for m in interesting_ligands if m.get("binding_quality")),
        key=lambda v: _tl_order.get(v, 3),
        default=None,
    )
    _grade, _lig_used = iridium_score(quality, _best_ligand_traffic)
    result["structure_quality"] = _grade
    result["structure_quality_ligand_used"] = _lig_used

    # Collect UniProt IDs: from input column first, then from resolved entities
    all_uniprot: list = []
    _input_uniprot = ""
    if uniprot_col:
        _raw = row.get(uniprot_col)
        _val = str(_raw).strip() if _raw is not None else ""
        if _val and _val.lower() != "nan":
            _input_uniprot = _val
            all_uniprot = [_val]
    for e in entities:
        for uid in e["uniprot_ids"]:
            if uid not in all_uniprot:
                all_uniprot.append(uid)

    if all_uniprot:
        if uniprot_col and not _input_uniprot:
            result[uniprot_col] = ",".join(all_uniprot)
            log.info("[%s] Resolved UniProt from PDB entry: %s", pdb_id, result[uniprot_col])
        elif not uniprot_col:
            result["resolved_uniprot"] = ",".join(all_uniprot)

    first_sequence = next((e["sequence"] for e in entities if e.get("sequence")), "")

    # Map each UniProt ID to the sequence length of its own entity so the
    # sibling/full-length split threshold is derived per-entity rather than
    # always from the first entity in the filtered list.
    _uid_seq_len: dict = {}
    for e in entities:
        if e.get("sequence"):
            for uid in e["uniprot_ids"]:
                if uid not in _uid_seq_len:
                    _uid_seq_len[uid] = len(e["sequence"])

    # --- Related entries ---
    fragment_ids: list = []
    sibling_ids: list = []
    fulllength_ids: list = []
    search_method = ""

    if all_uniprot:
        seen_fragments: set = set()
        seen_siblings: set = set()
        seen_fulllength: set = set()
        for uid in all_uniprot:
            uid_seq_len = _uid_seq_len.get(uid, len(first_sequence))
            if uid_seq_len:
                cache_key = (uid, uid_seq_len, _SEARCH_MAX_ROWS)
                if uniprot_search_cache is not None and cache_key in uniprot_search_cache:
                    log.info("[%s] UniProt search cache hit for %s", pdb_id, uid)
                    _frag, _sib, _full = uniprot_search_cache[cache_key]
                else:
                    log.info("[%s] Searching related entries by UniProt %s", pdb_id, uid)
                    _frag, _sib, _full = get_related_by_uniprot_split(
                        client, uid, uid_seq_len, max_rows=_SEARCH_MAX_ROWS
                    )
                    if uniprot_search_cache is not None:
                        uniprot_search_cache[cache_key] = (_frag, _sib, _full)
                for r in _frag:
                    if r not in seen_fragments:
                        seen_fragments.add(r)
                        fragment_ids.append(r)
                for r in _sib:
                    if r not in seen_siblings:
                        seen_siblings.add(r)
                        sibling_ids.append(r)
                for r in _full:
                    if r not in seen_fulllength:
                        seen_fulllength.add(r)
                        fulllength_ids.append(r)
                search_method = "uniprot_split"
            else:
                cache_key = (uid, 0, _SEARCH_MAX_ROWS)
                if uniprot_search_cache is not None and cache_key in uniprot_search_cache:
                    log.info("[%s] UniProt search cache hit for %s", pdb_id, uid)
                    _results = uniprot_search_cache[cache_key]
                else:
                    log.info("[%s] Searching related entries by UniProt %s", pdb_id, uid)
                    _results = get_related_by_uniprot(client, uid, max_rows=_SEARCH_MAX_ROWS)
                    if uniprot_search_cache is not None:
                        uniprot_search_cache[cache_key] = _results
                for r in _results:
                    if r not in seen_siblings:
                        seen_siblings.add(r)
                        sibling_ids.append(r)
                search_method = "uniprot"
    elif first_sequence:
        log.info("[%s] No UniProt; falling back to sequence similarity search", pdb_id)
        sibling_ids = get_related_by_sequence(
            client, first_sequence, identity_cutoff=seq_identity, max_rows=max_related
        )
        search_method = f"sequence_id_{seq_identity}"
    else:
        search_method = "none"

    # Remove self and other input CSV entries (their relationship is captured separately)
    other_input = input_pdb_ids - {pdb_id}
    fragment_ids = [r for r in fragment_ids if r.upper() not in other_input and r.upper() != pdb_id]
    sibling_ids = [r for r in sibling_ids if r.upper() not in other_input and r.upper() != pdb_id]
    fulllength_ids = [r for r in fulllength_ids if r.upper() not in other_input and r.upper() != pdb_id]

    def _get_related_data(pid: str) -> dict:
        if related_data_cache is not None:
            if pid not in related_data_cache:
                related_data_cache[pid] = _fetch_related_ligand_data(client, pid)
            return related_data_cache[pid]
        return _fetch_related_ligand_data(client, pid)

    # When include_all_related is set every found entry is fetched; otherwise cap at max_related.
    _fetch_limit = None if include_all_related else max_related

    fragment_no_ligand: list = []
    fragment_no_ligand_entries: list = []
    fragment_ligand_entries: list = []
    for frid in fragment_ids[:_fetch_limit]:
        log.info("[%s] Checking fragment %s for ligands", pdb_id, frid)
        data = _get_related_data(frid)
        if data["has_ligands"]:
            fragment_ligand_entries.append(data)
        else:
            fragment_no_ligand.append(frid)
            if include_all_related:
                fragment_no_ligand_entries.append(data)

    sibling_no_ligand: list = []
    sibling_no_ligand_entries: list = []
    sibling_ligand_entries: list = []
    for sid in sibling_ids[:_fetch_limit]:
        log.info("[%s] Checking sibling %s for ligands", pdb_id, sid)
        data = _get_related_data(sid)
        if data["has_ligands"]:
            sibling_ligand_entries.append(data)
        else:
            sibling_no_ligand.append(sid)
            if include_all_related:
                sibling_no_ligand_entries.append(data)

    fulllength_no_ligand: list = []
    fulllength_no_ligand_entries: list = []
    fulllength_ligand_entries: list = []
    for fid in fulllength_ids[:_fetch_limit]:
        log.info("[%s] Checking full-length %s for ligands", pdb_id, fid)
        data = _get_related_data(fid)
        if data["has_ligands"]:
            fulllength_ligand_entries.append(data)
        else:
            fulllength_no_ligand.append(fid)
            if include_all_related:
                fulllength_no_ligand_entries.append(data)

    result["all_fragment_pdb_ids"] = ",".join(fragment_ids)
    result["all_sibling_pdb_ids"] = ",".join(sibling_ids)
    result["all_fulllength_pdb_ids"] = ",".join(fulllength_ids)
    result["fragment_pdb_ids_no_ligand"] = ",".join(fragment_no_ligand)
    result["fragment_pdb_ids_no_ligand_count"] = len(fragment_no_ligand)
    result["sibling_pdb_ids_no_ligand"] = ",".join(sibling_no_ligand)
    result["sibling_pdb_ids_no_ligand_count"] = len(sibling_no_ligand)
    result["fulllength_pdb_ids_no_ligand"] = ",".join(fulllength_no_ligand)
    result["fulllength_pdb_ids_no_ligand_count"] = len(fulllength_no_ligand)
    result["related_search_method"] = search_method

    # Tag columns — None on primary rows; set on related-entry ligand sub-rows
    result["fragment_pdb_ids"] = None
    result["sibling_pdb_ids"] = None
    result["fulllength_pdb_ids"] = None

    # --- Binding sites ---
    all_site_features = []
    binding_sources: list = []

    for e in entities:
        all_site_features.extend(e.get("site_features") or [])

    pdbe_sites = get_pdbe_binding_sites(client, pdb_id) if valid_pdb else []

    # --- Known binders (still uses the first _MAX_RELATED_BINDER_ENTRIES siblings/full-lengths) ---
    query_ccd_codes = set(quality.get("ligands_present", "").split(",")) - {""}

    all_cofactors: list = []
    for e in entities:
        all_cofactors.extend(e.get("cofactors") or [])

    def _get_binders(pid: str) -> list:
        if binders_cache is not None:
            if pid not in binders_cache:
                binders_cache[pid] = extract_direct_binders(client, pid)
            return binders_cache[pid]
        return extract_direct_binders(client, pid)

    # Pull binders from fragments — unconditionally (fragments are subsets of the query domain)
    binder_fragments = [e["pdb_id"] for e in fragment_ligand_entries] + fragment_no_ligand
    for frid in binder_fragments[:_MAX_RELATED_BINDER_ENTRIES]:
        log.info("[%s] Fetching binders from fragment %s", pdb_id, frid)
        for b in _get_binders(frid):
            b = dict(b)
            b["binder_source_type"] = "fragment"
            all_cofactors.append(b)

    # Pull binders from siblings — prefer those with ligands first (already fetched), then no-ligand
    binder_siblings = [e["pdb_id"] for e in sibling_ligand_entries] + sibling_no_ligand
    for sid in binder_siblings[:_MAX_RELATED_BINDER_ENTRIES]:
        log.info("[%s] Fetching binders from sibling %s", pdb_id, sid)
        for b in _get_binders(sid):
            b = dict(b)
            b["binder_source_type"] = "sibling"
            all_cofactors.append(b)

    binder_fulllengths = [e["pdb_id"] for e in fulllength_ligand_entries] + fulllength_no_ligand
    for fid in binder_fulllengths[:_MAX_RELATED_BINDER_ENTRIES]:
        log.info("[%s] Fetching binders from full-length %s", pdb_id, fid)
        for b in _get_binders(fid):
            if b.get("chem_comp_id") and b["chem_comp_id"] in query_ccd_codes:
                b = dict(b)
                b["binder_source_type"] = "fulllength"
                all_cofactors.append(b)

    seen: set = set()
    unique_cofactors: list = []
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
        result["known_binders"] = ",".join(c["name"] for c in unique_cofactors if c.get("name"))
        result["known_binder_smiles"] = ",".join(
            c["smiles"] for c in unique_cofactors if c.get("smiles")
        )
    else:
        result["known_binders"] = ""
        result["known_binder_smiles"] = ""

    # --- Holo structure binding quality ---
    holo_results: list = []
    if unique_cofactors and all_uniprot:
        log.info("[%s] Looking up holo structures for %d known binder(s)", pdb_id, len(unique_cofactors))
        holo_results = get_holo_ligand_quality(client, unique_cofactors, all_uniprot)

    result["holo_quality"] = json.dumps(holo_results) if holo_results else ""

    # Primary row marker and flat ligand columns (None on primary rows)
    result["row_type"] = "primary"
    result["parent_pdb_id"] = None
    for col in LIGAND_DETAIL_COLS:
        result[col] = None

    # Internal metadata consumed by cli.py post-processing; stripped before CSV write
    result["_seq_len"] = len(first_sequence)
    result["_resolved_uniprot_ids"] = all_uniprot
    result["_ligand_metrics"] = ligand_metrics
    result["_peptide_entities"] = peptide_entities
    result["_fragment_ligand_entries"] = fragment_ligand_entries
    result["_fragment_no_ligand_entries"] = fragment_no_ligand_entries
    result["_sibling_ligand_entries"] = sibling_ligand_entries
    result["_sibling_no_ligand_entries"] = sibling_no_ligand_entries
    result["_fulllength_ligand_entries"] = fulllength_ligand_entries
    result["_fulllength_no_ligand_entries"] = fulllength_no_ligand_entries

    return result
