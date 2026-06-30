"""CLI entry point: parse arguments, iterate over input CSV, write enriched output."""

import argparse
import logging
import sys

import pandas as pd

from .client import RCSBClient
from .enrich import (
    LIGAND_DETAIL_COLS,
    _CRYSTAL_QUALITY_COLS,
    _normalise_pdb_id,
    build_ligand_rows,
    enrich_row,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

_PDB_COL_ALIASES = {"pdb_id", "pdbid", "pdb", "pdb id"}

_AUGMENTED_COLS = (
    # Row metadata
    "row_type",
    "parent_pdb_id",
    # Tag columns: single PDB ID set only on related-entry ligand sub-rows
    "related_pdb_ids",
    "fulllength_pdb_ids",
    "species",
    "structure_quality",
    "exp_method",
    # Structural quality (primary rows)
    "resolution_A",
    "r_work",
    "r_free",
    "clashscore",
    "ramachandran_outliers_pct",
    "rotamer_outliers_pct",
    "rsrz_outliers_pct",
    "bonds_rmsz",
    "angles_rmsz",
) + LIGAND_DETAIL_COLS + (  # flat ligand detail (ligand rows only)
    "ligands_present",
    "ligands_interesting",
    "ligands_noninteresting",
    "ligand_quality",
    "holo_quality",
    # Lists of related entries that carry NO meaningful ligand
    "related_pdb_ids_no_ligand",
    "related_pdb_ids_no_ligand_count",
    "fulllength_pdb_ids_no_ligand",
    "fulllength_pdb_ids_no_ligand_count",
    "related_search_method",
    "binding_site_sources",
    "binding_site_notes",
    "known_binders",
    "known_binder_smiles",
)

_AUGMENTED_COLS_SET = frozenset(_AUGMENTED_COLS)

# Internal metadata keys stripped before writing CSV
_INTERNAL_KEYS = {
    "_seq_len", "_resolved_uniprot_ids",
    "_ligand_metrics", "_peptide_entities",
    "_sibling_ligand_entries", "_fulllength_ligand_entries",
}


def detect_pdb_col(columns: list, hint: str | None) -> str:
    if hint:
        if hint in columns:
            return hint
        raise ValueError(f"Column '{hint}' not found. Available: {columns}")
    for col in columns:
        if col.strip().lower() in _PDB_COL_ALIASES:
            return col
    raise ValueError(f"Cannot auto-detect PDB ID column. Use --pdb-col. Available: {columns}")


def detect_uniprot_col(columns: list, hint: str | None) -> str | None:
    if hint:
        if hint in columns:
            return hint
        log.warning("UniProt column '%s' not found; skipping UniProt lookup", hint)
        return None
    for col in columns:
        if col.strip().lower() in {"uniprot", "uniprot_id", "uniprot_acc", "accession"}:
            return col
    return None


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
    input_cols = df.columns.tolist()

    log.info(
        "Input: %d rows | PDB column: '%s' | UniProt column: %s",
        len(df),
        pdb_col,
        uniprot_col or "none",
    )

    input_pdb_ids = frozenset(
        _normalise_pdb_id(str(v).strip()).upper()
        for v in df[pdb_col]
        if str(v).strip()
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
            log.warning(
                "Row %d: normalised PDB ID '%s' → '%s' (Excel formatting?)",
                i, row.get(pdb_col, ""), pdb_id,
            )
        log.info("Processing row %d/%d: %s", i, len(df), pdb_id)
        try:
            enriched = enrich_row(
                row=row,
                client=client,
                pdb_col=pdb_col,
                uniprot_col=uniprot_col,
                seq_identity=args.seq_identity,
                max_related=args.max_related,
                input_pdb_ids=input_pdb_ids,
            )
        except Exception as exc:
            log.error("[%s] Failed: %s", pdb_id, exc, exc_info=True)
            enriched = dict(row)
            for col in _AUGMENTED_COLS_SET:
                enriched.setdefault(col, None)
            enriched["_ligand_metrics"] = []
            enriched["_peptide_entities"] = []
            enriched["_sibling_ligand_entries"] = []
            enriched["_fulllength_ligand_entries"] = []
            enriched["_seq_len"] = 0
            enriched["_resolved_uniprot_ids"] = []
        enriched["_pdb_id"] = pdb_id
        enriched_rows.append(enriched)

    all_output_cols = tuple(input_cols) + _AUGMENTED_COLS
    final_rows = []
    # Track which (parent_pdb_id, source_pdb_id, ligand_id, chain_id) tuples have been emitted
    # to suppress exact duplicate sub-rows from duplicate input PDB IDs or repeated chains.
    emitted_ligand_keys: set = set()

    def _add_ligand_rows(sub_rows: list) -> None:
        for r in sub_rows:
            key = (
                r.get("parent_pdb_id"),
                r.get("related_pdb_ids"),
                r.get("fulllength_pdb_ids"),
                r.get("ligand_id"),
                r.get("chain_id"),
            )
            if key not in emitted_ligand_keys:
                emitted_ligand_keys.add(key)
                final_rows.append(r)

    for row in enriched_rows:
        pdb_id_for_ligands = row.pop("_pdb_id", "")
        sibling_entries = row.pop("_sibling_ligand_entries", []) or []
        fulllength_entries = row.pop("_fulllength_ligand_entries", []) or []
        for key in _INTERNAL_KEYS:
            row.pop(key, None)

        for col in _AUGMENTED_COLS:
            row.setdefault(col, None)
        final_rows.append(row)

        # Sibling ligands — tagged with the sibling's PDB ID and its crystal quality
        for entry in sibling_entries:
            _add_ligand_rows(build_ligand_rows(
                pdb_id=pdb_id_for_ligands,
                ligand_metrics=entry["ligand_metrics"],
                peptide_entities=entry["peptide_entities"],
                all_output_cols=all_output_cols,
                tags={"related_pdb_ids": entry["pdb_id"], "species": entry.get("species", ""), "structure_quality": entry.get("structure_quality", ""), **entry.get("entry_quality", {})},
            ))

        # Full-length ligands — tagged in fulllength_pdb_ids with its crystal quality
        for entry in fulllength_entries:
            _add_ligand_rows(build_ligand_rows(
                pdb_id=pdb_id_for_ligands,
                ligand_metrics=entry["ligand_metrics"],
                peptide_entities=entry["peptide_entities"],
                all_output_cols=all_output_cols,
                tags={"fulllength_pdb_ids": entry["pdb_id"], "species": entry.get("species", ""), "structure_quality": entry.get("structure_quality", ""), **entry.get("entry_quality", {})},
            ))

    out_df = pd.DataFrame(final_rows, columns=list(all_output_cols))
    out_df.to_csv(args.output, index=False)
    log.info("Written %d rows (%d primary) to %s", len(out_df), len(enriched_rows), args.output)
