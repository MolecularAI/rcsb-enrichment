"""Related PDB entry search: by UniProt (with optional length split) or sequence similarity."""

from .client import RCSBClient, RCSB_SEARCH

# Chains longer than query × this factor are "full-length" proteins.
_SEQ_LEN_RATIO = 1.4


def _uniprot_length_search(
    client: RCSBClient,
    uniprot_id: str,
    length_operator: str,
    length_value: int,
    max_rows: int,
) -> list:
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
    """Split same-UniProt entries into (siblings, full_length) by chain length.

    'Siblings' have chain length ≤ query × _SEQ_LEN_RATIO.
    'Full-length' have chain length > that threshold.
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
