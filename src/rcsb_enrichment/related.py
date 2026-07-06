"""Related PDB entry search: by UniProt (with optional length split) or sequence similarity."""

from .client import RCSBClient, RCSB_SEARCH

# Chains longer than query × this factor are "full-length" proteins.
_SEQ_LEN_RATIO_UPPER = 1.4
# Chains shorter than query × this factor are "fragment" (subset) proteins.
_SEQ_LEN_RATIO_LOWER = 0.8


def _uniprot_length_search(
    client: RCSBClient,
    uniprot_id: str,
    length_operator: str,
    length_value: int,
    max_rows: int,
) -> list:
    # Use rcsb_entry_info.deposited_polymer_monomer_count (total residues across ALL protein
    # chains in the entry) as the size signal.  A per-chain search on
    # entity_poly.rcsb_sample_sequence_length would miss structural subsets like 1JFF vs 5S5V:
    # both have the same tubulin chain length (451 aa), but 1JFF has only 2 protein chains while
    # 5S5V has 4, so 1JFF's total is ~896 vs 5S5V's ~2319 — the entry-level count captures this.
    # return_type="entry" is correct here because deposited_polymer_monomer_count is an
    # entry-level attribute; the UniProt condition selects only entries that contain at least
    # one chain matching the accession.
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
                        "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
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
    return [r["identifier"] for r in (data or {}).get("result_set") or []]


def get_related_by_uniprot(client: RCSBClient, uniprot_id: str, max_rows: int = 200) -> list:
    """All entries sharing uniprot_id (no length split — used when query seq len is unknown)."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers"
                             ".reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows}},
    }
    data = client.post(RCSB_SEARCH, payload)
    return [r["identifier"] for r in (data or {}).get("result_set") or []]


def get_related_by_uniprot_split(
    client: RCSBClient,
    uniprot_id: str,
    query_seq_len: int,
    max_rows: int = 200,
) -> tuple:
    """Split same-UniProt entries into (fragments, siblings, full_length) by chain length.

    'Fragments' have chain length < query × _SEQ_LEN_RATIO_LOWER  (subsets of the query).
    'Siblings'  have chain length in [lower_threshold, upper_threshold].
    'Full-length' have chain length > query × _SEQ_LEN_RATIO_UPPER.
    Three search requests; no per-hit extra fetches.
    Returns (fragments, siblings, full_length).
    """
    lower = int(query_seq_len * _SEQ_LEN_RATIO_LOWER)
    upper = int(query_seq_len * _SEQ_LEN_RATIO_UPPER)
    fragments = _uniprot_length_search(client, uniprot_id, "less", lower, max_rows)
    # Siblings: entries with length in [lower, upper] — use the full ≤upper set minus fragments.
    # A single range search is not available in the RCSB text service, so we fetch ≤upper and
    # subtract the fragment set, which avoids a fourth request.
    lte_upper = _uniprot_length_search(client, uniprot_id, "less_or_equal", upper, max_rows)
    fragment_set = set(fragments)
    siblings = [r for r in lte_upper if r not in fragment_set]
    full_length = _uniprot_length_search(client, uniprot_id, "greater", upper, max_rows)
    return fragments, siblings, full_length


def get_related_by_sequence(
    client: RCSBClient,
    sequence: str,
    identity_cutoff: float = 0.9,
    max_rows: int = 25,
) -> list:
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
    return [r["identifier"] for r in (data or {}).get("result_set") or []]
