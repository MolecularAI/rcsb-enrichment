"""Tests for holo._inchikey_to_ccd, _find_holo_entries, get_holo_ligand_quality."""

import pytest
from unittest.mock import MagicMock, call
from rcsb_enrichment.holo import (
    _MAX_HOLO_ENTRIES,
    _find_holo_entries,
    _inchikey_to_ccd,
    get_holo_ligand_quality,
)


def _search_result(*identifiers):
    return {"result_set": [{"identifier": i} for i in identifiers]}


# ---------------------------------------------------------------------------
# _inchikey_to_ccd
# ---------------------------------------------------------------------------

class TestInchiKeyToCcd:
    def test_resolves_inchikey(self):
        client = MagicMock()
        client.post.return_value = _search_result("ATP")
        assert _inchikey_to_ccd(client, "ZKHQWZAMYRWXGA-KQYNXXCUSA-N") == "ATP"

    def test_no_result_returns_none(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        assert _inchikey_to_ccd(client, "NOTANINCHIKEY") is None

    def test_none_response_returns_none(self):
        client = MagicMock()
        client.post.return_value = None
        assert _inchikey_to_ccd(client, "SOMEKEY") is None

    def test_empty_inchikey_skips_api(self):
        client = MagicMock()
        assert _inchikey_to_ccd(client, "") is None
        assert _inchikey_to_ccd(client, None) is None
        client.post.assert_not_called()

    def test_query_uses_text_chem_service(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _inchikey_to_ccd(client, "SOMEKEY")
        payload = client.post.call_args[0][1]
        assert payload["query"]["service"] == "text_chem"

    def test_query_uses_exact_match(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _inchikey_to_ccd(client, "SOMEKEY")
        payload = client.post.call_args[0][1]
        assert payload["query"]["parameters"]["operator"] == "exact_match"
        assert payload["query"]["parameters"]["value"] == "SOMEKEY"

    def test_return_type_is_mol_definition(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _inchikey_to_ccd(client, "SOMEKEY")
        payload = client.post.call_args[0][1]
        assert payload["return_type"] == "mol_definition"


# ---------------------------------------------------------------------------
# _find_holo_entries
# ---------------------------------------------------------------------------

class TestFindHoloEntries:
    def test_returns_pdb_ids(self):
        client = MagicMock()
        client.post.return_value = _search_result("1ABC", "2DEF")
        result = _find_holo_entries(client, "P12345", "ATP")
        assert result == ["1ABC", "2DEF"]

    def test_no_results(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        assert _find_holo_entries(client, "P12345", "ATP") == []

    def test_none_response(self):
        client = MagicMock()
        client.post.return_value = None
        assert _find_holo_entries(client, "P12345", "ATP") == []

    def test_query_includes_uniprot_and_ccd(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _find_holo_entries(client, "P12345", "LIG")
        payload = client.post.call_args[0][1]
        payload_str = str(payload)
        assert "P12345" in payload_str
        assert "LIG" in payload_str

    def test_return_type_is_entry(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _find_holo_entries(client, "P12345", "LIG")
        payload = client.post.call_args[0][1]
        assert payload["return_type"] == "entry"

    def test_respects_max_rows(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _find_holo_entries(client, "P12345", "ATP", max_rows=2)
        payload = client.post.call_args[0][1]
        assert payload["request_options"]["paginate"]["rows"] == 2

    def test_ccd_attribute_is_correct(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        _find_holo_entries(client, "P12345", "ATP")
        payload = client.post.call_args[0][1]
        nodes = payload["query"]["nodes"]
        ccd_node = next(n for n in nodes if "nonpolymer_comp_id" in str(n))
        assert ccd_node["parameters"]["attribute"] == (
            "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"
        )


# ---------------------------------------------------------------------------
# get_holo_ligand_quality
# ---------------------------------------------------------------------------

class TestGetHoloLigandQuality:
    def _entry_data(self, np_ids):
        return {
            "rcsb_entry_container_identifiers": {
                "non_polymer_entity_ids": np_ids
            }
        }

    def _entity_data(self, ccd, asym_ids):
        return {
            "rcsb_nonpolymer_entity_container_identifiers": {
                "nonpolymer_comp_id": ccd,
                "asym_ids": asym_ids,
            },
            "rcsb_nonpolymer_entity": {"pdbx_description": ccd},
        }

    def _inst_data(self, rscc=0.85, isi=True):
        return {
            "rcsb_nonpolymer_instance_validation_score": [{
                "RSCC": rscc, "RSR": 0.15,
                "mogul_bonds_RMSZ": 1.0, "mogul_angles_RMSZ": 1.0,
                "intermolecular_clashes": 0,
                "is_subject_of_investigation": isi,
            }],
            "rcsb_target_neighbors": [],
        }

    def test_no_uniprot_ids_returns_empty(self):
        client = MagicMock()
        assert get_holo_ligand_quality(client, [{"name": "Drug", "chem_comp_id": "LIG"}], []) == []

    def test_binder_with_direct_ccd(self):
        client = MagicMock()
        # post: find_holo_entries → returns 1ABC
        client.post.return_value = _search_result("1ABC")
        # get: entry → entity → instance
        client.get.side_effect = [
            self._entry_data(["1"]),
            self._entity_data("LIG", ["A"]),
            self._inst_data(),
        ]
        results = get_holo_ligand_quality(
            client,
            [{"name": "TestDrug", "chem_comp_id": "LIG", "inchikey": "SOMEKEY"}],
            ["P12345"],
        )
        assert len(results) == 1
        assert results[0]["ccd_code"] == "LIG"
        assert results[0]["binder_name"] == "TestDrug"
        assert "1ABC" in results[0]["holo_pdb_ids"]
        assert len(results[0]["holo_ligand_metrics"]) == 1

    def test_binder_inchikey_resolved_to_ccd(self):
        client = MagicMock()
        # first post: _inchikey_to_ccd → "GTP"; second post: _find_holo_entries → "2GTP"
        client.post.side_effect = [
            _search_result("GTP"),    # InChIKey resolution
            _search_result("2GTP"),   # holo search
        ]
        client.get.side_effect = [
            self._entry_data(["1"]),
            self._entity_data("GTP", ["A"]),
            self._inst_data(),
        ]
        results = get_holo_ligand_quality(
            client,
            [{"name": "GTP", "chem_comp_id": None, "inchikey": "ZKHQWZAMYRWXGA-KQYNXXCUSA-N"}],
            ["P12345"],
        )
        assert len(results) == 1
        assert results[0]["ccd_code"] == "GTP"

    def test_no_ccd_resolvable_skipped(self):
        client = MagicMock()
        # InChIKey resolution returns nothing
        client.post.return_value = {"result_set": []}
        results = get_holo_ligand_quality(
            client,
            [{"name": "Unknown", "chem_comp_id": None, "inchikey": ""}],
            ["P12345"],
        )
        assert results == []

    def test_no_holo_entries_found_skipped(self):
        client = MagicMock()
        # CCD known, but holo search returns no entries
        client.post.return_value = {"result_set": []}
        results = get_holo_ligand_quality(
            client,
            [{"name": "Drug", "chem_comp_id": "LIG", "inchikey": ""}],
            ["P12345"],
        )
        assert results == []

    def test_holo_entry_404_skipped(self):
        client = MagicMock()
        client.post.return_value = _search_result("1ABC")
        client.get.return_value = None  # 404 on entry fetch
        results = get_holo_ligand_quality(
            client,
            [{"name": "Drug", "chem_comp_id": "LIG", "inchikey": ""}],
            ["P12345"],
        )
        assert results == []

    def test_only_matching_ccd_ligands_included(self):
        client = MagicMock()
        client.post.return_value = _search_result("1ABC")
        # Entry has two nonpolymer entities: LIG and HOH
        client.get.side_effect = [
            self._entry_data(["1", "2"]),
            self._entity_data("LIG", ["A"]),
            self._inst_data(rscc=0.88),
            self._entity_data("HOH", ["W"]),
            self._inst_data(rscc=0.50, isi=False),
        ]
        results = get_holo_ligand_quality(
            client,
            [{"name": "Drug", "chem_comp_id": "LIG", "inchikey": ""}],
            ["P12345"],
        )
        assert len(results) == 1
        assert all(m["ligand_id"] == "LIG" for m in results[0]["holo_ligand_metrics"])

    def test_multiple_binders_each_resolved(self):
        client = MagicMock()
        # Two binders with direct CCD codes; each search finds one entry
        client.post.side_effect = [
            _search_result("1AB1"),  # holo for LIG1
            _search_result("1AB2"),  # holo for LIG2
        ]
        def _get(url):
            if "1AB1" in url:
                return self._entry_data(["1"]) if url.endswith("1AB1") else self._entity_data("LIG1", ["A"])
            if "1AB2" in url:
                return self._entry_data(["1"]) if url.endswith("1AB2") else self._entity_data("LIG2", ["A"])
            return self._inst_data()
        client.get.side_effect = [
            self._entry_data(["1"]),    # entry for 1AB1
            self._entity_data("LIG1", ["A"]),
            self._inst_data(),
            self._entry_data(["1"]),    # entry for 1AB2
            self._entity_data("LIG2", ["A"]),
            self._inst_data(),
        ]
        binders = [
            {"name": "Drug1", "chem_comp_id": "LIG1", "inchikey": ""},
            {"name": "Drug2", "chem_comp_id": "LIG2", "inchikey": ""},
        ]
        results = get_holo_ligand_quality(client, binders, ["P12345"])
        assert len(results) == 2
        ccd_codes = {r["ccd_code"] for r in results}
        assert ccd_codes == {"LIG1", "LIG2"}

    def test_uses_first_uniprot_id(self):
        client = MagicMock()
        client.post.return_value = {"result_set": []}
        get_holo_ligand_quality(
            client,
            [{"name": "Drug", "chem_comp_id": "LIG", "inchikey": ""}],
            ["P11111", "P22222"],
        )
        payload = client.post.call_args[0][1]
        assert "P11111" in str(payload)
        assert "P22222" not in str(payload)

    def test_max_holo_entries_cap_respected(self):
        client = MagicMock()
        # holo search returns more entries than _MAX_HOLO_ENTRIES
        many_ids = [f"{i:04d}" for i in range(1, _MAX_HOLO_ENTRIES + 5)]
        client.post.return_value = _search_result(*many_ids)
        client.get.return_value = None  # all 404 — just test call count
        get_holo_ligand_quality(
            client,
            [{"name": "Drug", "chem_comp_id": "LIG", "inchikey": ""}],
            ["P12345"],
        )
        # find_holo_entries itself passes max_rows=_MAX_HOLO_ENTRIES to the search,
        # so the search will already return at most _MAX_HOLO_ENTRIES results
        payload = client.post.call_args[0][1]
        assert payload["request_options"]["paginate"]["rows"] == _MAX_HOLO_ENTRIES
