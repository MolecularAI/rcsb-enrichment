"""Tests for related.get_related_by_uniprot*, get_related_by_sequence."""

import pytest
from unittest.mock import MagicMock, call
from rcsb_enrichment.related import (
    _SEQ_LEN_RATIO_LOWER,
    _SEQ_LEN_RATIO_UPPER,
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
    def test_splits_into_fragments_siblings_full_length(self):
        client = MagicMock()
        client.post.side_effect = [
            _search_result("FRAG1"),             # fragments (less)
            _search_result("SIB1", "SIB2"),      # lte_upper set (less_or_equal)
            _search_result("FULL1"),             # full-length (greater)
        ]
        fragments, siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert fragments == ["FRAG1"]
        assert siblings == ["SIB1", "SIB2"]
        assert full == ["FULL1"]

    def test_fragments_excluded_from_siblings(self):
        client = MagicMock()
        # FRAG1 appears in both the less search and the less_or_equal search;
        # it must be subtracted so siblings only contains the non-fragment entries.
        client.post.side_effect = [
            _search_result("FRAG1"),             # fragments
            _search_result("FRAG1", "SIB1"),     # lte_upper — includes the fragment
            _search_result(),                    # full-length
        ]
        fragments, siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert "FRAG1" not in siblings
        assert "SIB1" in siblings

    def test_upper_threshold_is_1_4x(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_uniprot_split(client, "P12345", 200)
        calls = client.post.call_args_list
        # call[0]: less (lower threshold = int(200 * 0.8) = 160)
        # call[1]: less_or_equal (upper threshold = int(200 * 1.4) = 280)
        # call[2]: greater (upper threshold = int(200 * 1.4) = 280)
        upper = int(200 * _SEQ_LEN_RATIO_UPPER)
        lower = int(200 * _SEQ_LEN_RATIO_LOWER)

        first_payload = calls[0][0][1]
        assert first_payload["query"]["nodes"][1]["parameters"]["operator"] == "less"
        assert first_payload["query"]["nodes"][1]["parameters"]["value"] == lower

        second_payload = calls[1][0][1]
        assert second_payload["query"]["nodes"][1]["parameters"]["operator"] == "less_or_equal"
        assert second_payload["query"]["nodes"][1]["parameters"]["value"] == upper

        third_payload = calls[2][0][1]
        assert third_payload["query"]["nodes"][1]["parameters"]["operator"] == "greater"
        assert third_payload["query"]["nodes"][1]["parameters"]["value"] == upper

    def test_three_api_calls_made(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_related_by_uniprot_split(client, "P12345", 100)
        assert client.post.call_count == 3

    def test_all_empty(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        fragments, siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert fragments == []
        assert siblings == []
        assert full == []

    def test_returns_three_tuple(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        result = get_related_by_uniprot_split(client, "P12345", 100)
        assert len(result) == 3

    def test_none_response_handled(self):
        client = MagicMock()
        client.post.return_value = None
        fragments, siblings, full = get_related_by_uniprot_split(client, "P12345", 100)
        assert fragments == []
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
