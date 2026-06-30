"""Tests for related.get_related_by_uniprot*, get_related_by_sequence."""

import pytest
from unittest.mock import MagicMock, call
from rcsb_enrichment.related import (
    _SEQ_LEN_RATIO,
    get_related_by_sequence,
    get_related_by_uniprot,
    get_related_by_uniprot_split,
)


def _search_result(*identifiers):
    return {"result_set": [{"identifier": i} for i in identifiers]}


class TestGetRelatedByUniprot:
    def test_returns_identifiers(self):
        client = MagicMock()
        client.post.return_value = _search_result("1ABC", "2DEF")
        result = get_related_by_uniprot(client, "P12345")
        assert result == ["1ABC", "2DEF"]

    def test_empty_result(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        assert get_related_by_uniprot(client, "P12345") == []

    def test_none_response(self):
        client = MagicMock()
        client.post.return_value = None
        assert get_related_by_uniprot(client, "P12345") == []

    def test_query_uses_exact_match(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_uniprot(client, "P12345")
        payload = client.post.call_args[0][1]
        assert payload["query"]["parameters"]["operator"] == "exact_match"
        assert payload["query"]["parameters"]["value"] == "P12345"
        assert payload["return_type"] == "entry"


class TestGetRelatedByUniprotSplit:
    def test_splits_into_siblings_and_full_length(self):
        client = MagicMock()
        client.post.side_effect = [
            _search_result("SIB1", "SIB2"),   # siblings (less_or_equal)
            _search_result("FULL1"),            # full-length (greater)
        ]
        siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert siblings == ["SIB1", "SIB2"]
        assert full == ["FULL1"]

    def test_threshold_is_1_4x(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_uniprot_split(client, "P12345", 200)
        calls = client.post.call_args_list
        # First call: less_or_equal, threshold = int(200 * 1.4) = 280
        first_payload = calls[0][0][1]
        length_node = first_payload["query"]["nodes"][1]
        assert length_node["parameters"]["operator"] == "less_or_equal"
        assert length_node["parameters"]["value"] == int(200 * _SEQ_LEN_RATIO)
        # Second call: greater, same threshold
        second_payload = calls[1][0][1]
        length_node2 = second_payload["query"]["nodes"][1]
        assert length_node2["parameters"]["operator"] == "greater"

    def test_two_api_calls_made(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_uniprot_split(client, "P12345", 100)
        assert client.post.call_count == 2

    def test_both_empty(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert siblings == []
        assert full == []


class TestGetRelatedBySequence:
    def test_returns_identifiers(self):
        client = MagicMock()
        client.post.return_value = _search_result("1ABC", "2DEF")
        result = get_related_by_sequence(client, "ACDEFGHIK")
        assert result == ["1ABC", "2DEF"]

    def test_uses_correct_identity_cutoff(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_sequence(client, "ACDEFGHIK", identity_cutoff=0.95)
        payload = client.post.call_args[0][1]
        assert payload["query"]["parameters"]["identity_cutoff"] == 0.95
        assert payload["query"]["parameters"]["sequence_type"] == "protein"

    def test_none_response(self):
        client = MagicMock()
        client.post.return_value = None
        assert get_related_by_sequence(client, "ACDEFGHIK") == []
