"""Tests for enrich._normalise_pdb_id, _is_valid_pdb_id, build_ligand_rows, and enrich_row."""

import json
import pytest
from unittest.mock import MagicMock, patch
from rcsb_enrichment.enrich import (
    LIGAND_DETAIL_COLS,
    _collect_entity_names,
    _entity_matches_names,
    _fetch_related_ligand_data,
    _is_valid_pdb_id,
    _normalise_pdb_id,
    build_ligand_rows,
    enrich_row,
)
from rcsb_enrichment.cli import _AUGMENTED_COLS


# ---------------------------------------------------------------------------
# _normalise_pdb_id
# ---------------------------------------------------------------------------

class TestNormalisePdbId:
    def test_normal_id_unchanged(self):
        assert _normalise_pdb_id("1ABC") == "1ABC"

    def test_lowercase_preserved(self):
        # normalise doesn't uppercase; caller does
        assert _normalise_pdb_id("1abc") == "1abc"

    def test_excel_thousands_separator(self):
        # Excel formats "6BHD" as "6,000 BHD" → strip commas and spaces → "6000BHD" → take 4 → "6000"
        # Actual mangling: "6,BHD" is more likely for short IDs; general rule: strip , and spaces
        assert _normalise_pdb_id("6,BHD") == "6BHD"

    def test_spaces_stripped(self):
        assert _normalise_pdb_id(" 1 ABC ") == "1ABC"

    def test_commas_stripped(self):
        assert _normalise_pdb_id("1,ABC") == "1ABC"

    def test_too_short_returned_as_is(self):
        assert _normalise_pdb_id("AB") == "AB"

    def test_truncated_to_four(self):
        assert _normalise_pdb_id("1ABCDEF") == "1ABC"


class TestIsValidPdbId:
    def test_valid_ids(self):
        for pid in ("1ABC", "4Z02", "6BHD", "1000", "AAAA"):
            assert _is_valid_pdb_id(pid), f"{pid} should be valid"

    def test_lowercase_valid(self):
        assert _is_valid_pdb_id("1abc")

    def test_too_short(self):
        assert not _is_valid_pdb_id("1AB")

    def test_too_long(self):
        assert not _is_valid_pdb_id("1ABCD")

    def test_invalid_chars(self):
        assert not _is_valid_pdb_id("1AB!")
        assert not _is_valid_pdb_id("1AB-")

    def test_alphafold_passes_regex(self):
        # "ALPH" is technically 4 alphanumeric chars — passes the regex.
        # The guard against non-PDB IDs is downstream: get_entry_quality returns {} on 404.
        assert _is_valid_pdb_id("ALPH")


# ---------------------------------------------------------------------------
# build_ligand_rows
# ---------------------------------------------------------------------------

class TestBuildLigandRows:
    def _all_cols(self):
        return ("pdb_id", "uniprot") + _AUGMENTED_COLS

    def _good_metric(self, ligand_id="LIG", chain="A"):
        return {
            "ligand_id": ligand_id,
            "chain_id": chain,
            "is_interesting": True,
            "contact_residue_count": 5,
            "binding_quality": "good",
            "ligand_rscc": 0.90,
            "ligand_rsr": 0.12,
            "ligand_rmsz_bonds": 1.0,
            "ligand_rmsz_angles": 1.0,
            "ligand_intermolecular_clashes": 0,
            "contact_outlier_fraction": 0.05,
            "contact_residues": "ALA10(A)",
        }

    def test_interesting_ligand_with_contacts_included(self):
        rows = build_ligand_rows("1ABC", [self._good_metric()], [], self._all_cols())
        assert len(rows) == 1
        assert rows[0]["row_type"] == "ligand"
        assert rows[0]["ligand_type"] == "small_molecule"
        assert rows[0]["ligand_id"] == "LIG"
        assert rows[0]["parent_pdb_id"] == "1ABC"

    def test_non_interesting_excluded(self):
        m = self._good_metric()
        m["is_interesting"] = False
        assert build_ligand_rows("1ABC", [m], [], self._all_cols()) == []

    def test_interesting_no_contacts_no_quality_excluded(self):
        m = self._good_metric()
        m["contact_residue_count"] = 0
        m["binding_quality"] = ""
        assert build_ligand_rows("1ABC", [m], [], self._all_cols()) == []

    def test_interesting_no_contacts_but_has_quality_included(self):
        m = self._good_metric()
        m["contact_residue_count"] = 0
        m["binding_quality"] = "fair"
        rows = build_ligand_rows("1ABC", [m], [], self._all_cols())
        assert len(rows) == 1

    def test_peptide_always_included(self):
        peptides = [{"bird_id": "PRD_000001", "sequence": "ACDEF"}]
        rows = build_ligand_rows("1ABC", [], peptides, self._all_cols())
        assert len(rows) == 1
        assert rows[0]["ligand_type"] == "peptide"
        assert rows[0]["ligand_id"] == "PRD_000001"

    def test_peptide_sequence_used_when_no_bird_id(self):
        peptides = [{"bird_id": None, "sequence": "ACDEFGH"}]
        rows = build_ligand_rows("1ABC", [], peptides, self._all_cols())
        assert rows[0]["ligand_id"] == "ACDEFGH"

    def test_empty_peptide_label_skipped(self):
        peptides = [{"bird_id": None, "sequence": ""}]
        assert build_ligand_rows("1ABC", [], peptides, self._all_cols()) == []

    def test_input_cols_are_none_on_ligand_rows(self):
        rows = build_ligand_rows("1ABC", [self._good_metric()], [], self._all_cols())
        # 'pdb_id' and 'uniprot' are input passthrough cols → must be None
        assert rows[0]["pdb_id"] is None
        assert rows[0]["uniprot"] is None

    def test_multiple_ligands(self):
        metrics = [self._good_metric("LIG1", "A"), self._good_metric("LIG2", "B")]
        rows = build_ligand_rows("1ABC", metrics, [], self._all_cols())
        assert len(rows) == 2
        assert {r["ligand_id"] for r in rows} == {"LIG1", "LIG2"}


# ---------------------------------------------------------------------------
# enrich_row (integration-style, all API calls mocked)
# ---------------------------------------------------------------------------

def _make_client(*, quality=None, entities=None, ligands=None, related=None):
    """Build a mock client with pre-configured side_effect responses."""
    client = MagicMock()
    client.get.return_value = None
    client.post.return_value = {"result_set": []}
    return client


class TestEnrichRow:
    def _base_row(self, pdb_id="1ABC", uniprot="P12345"):
        return {"PDBID": pdb_id, "Uniprot": uniprot}

    def _mock_full_run(self, pdb_id="1ABC"):
        """Return a client mock that satisfies a minimal full enrich_row call."""
        entry_data = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {
                "diffrn_resolution_high": {"value": 2.0},
                "nonpolymer_bound_components": [],
            },
            "refine": [{"ls_R_factor_R_work": 0.18, "ls_R_factor_R_free": 0.22}],
            "pdbx_vrpt_summary_geometry": [{"clashscore": 1.5, "percent_ramachandran_outliers": 0.0,
                                             "percent_rotamer_outliers": 0.0, "bonds_RMSZ": 0.9, "angles_RMSZ": 1.0}],
            "pdbx_vrpt_summary_diffraction": [{"percent_RSRZ_outliers": 2.0}],
            "rcsb_entry_container_identifiers": {
                "polymer_entity_ids": ["1"],
                "non_polymer_entity_ids": [],
            },
        }
        entity_data = {
            "entity_poly": {"rcsb_entity_polymer_type": "Protein", "pdbx_seq_one_letter_code_can": "A" * 100},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P12345"], "bird_id": None},
            "rcsb_target_cofactors": [],
            "rcsb_polymer_entity_feature": [],
        }
        pdbe_data = {pdb_id.lower(): []}
        client = MagicMock()
        client.get.side_effect = [entry_data, entity_data, pdbe_data]
        client.post.return_value = {"result_set": []}
        return client

    def test_basic_fields_populated(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert result["exp_method"] == "X-RAY DIFFRACTION"
        assert result["resolution_A"] == pytest.approx(2.0)
        assert result["row_type"] == "primary"
        assert result["parent_pdb_id"] is None

    def test_internal_metadata_attached(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert "_seq_len" in result
        assert "_resolved_uniprot_ids" in result
        assert "_ligand_metrics" in result
        assert "_peptide_entities" in result

    def test_seq_len_matches_sequence(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert result["_seq_len"] == 100

    def test_input_pdb_excluded_from_related(self):
        client = self._mock_full_run()
        # get() side_effect: entry, entity, pdbe — then None for each _fetch_related_ligand_data call
        client.get.side_effect = list(client.get.side_effect) + [None, None]
        client.post.side_effect = [
            {"result_set": [{"identifier": "1ABC"}, {"identifier": "2DEF"}]},  # siblings
            {"result_set": []},                                                   # full-length
        ]
        with patch("rcsb_enrichment.enrich.extract_direct_binders", return_value=[]):
            result = enrich_row(
                row=self._base_row(),
                client=client,
                pdb_col="PDBID",
                uniprot_col="Uniprot",
                seq_identity=0.9,
                max_related=25,
                input_pdb_ids=frozenset(["1ABC"]),
            )
        assert "1ABC" not in result.get("related_pdb_ids_no_ligand", "")

    def test_other_input_pdb_excluded_from_related(self):
        client = self._mock_full_run()
        # get() side_effect: entry, entity, pdbe — plus one for _fetch_related_ligand_data(3GHI)
        client.get.side_effect = list(client.get.side_effect) + [None]
        client.post.side_effect = [
            {"result_set": [{"identifier": "2DEF"}, {"identifier": "3GHI"}]},  # siblings
            {"result_set": []},
        ]
        with patch("rcsb_enrichment.enrich.extract_direct_binders", return_value=[]):
            result = enrich_row(
                row=self._base_row(),
                client=client,
                pdb_col="PDBID",
                uniprot_col="Uniprot",
                seq_identity=0.9,
                max_related=25,
                input_pdb_ids=frozenset(["1ABC", "2DEF"]),
            )
        assert "2DEF" not in result.get("related_pdb_ids_no_ligand", "")
        assert "3GHI" in result.get("related_pdb_ids_no_ligand", "")

    def test_404_entry_gives_none_quality_fields(self):
        # When get() returns None (404) quality = {} → all quality fields are None
        client = MagicMock()
        client.get.return_value = None
        client.post.return_value = {"result_set": []}
        result = enrich_row(
            row={"PDBID": "1ABC", "Uniprot": "P12345"},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert result["exp_method"] is None
        assert result["resolution_A"] is None

    def test_all_resolved_uniprots_searched_for_related_entries(self):
        """Bug: when CSV has no UniProt and the entry has multiple protein chains with
        different UniProt IDs, the related-entry search must be performed for ALL resolved
        IDs, not just all_uniprot[0].  Previously only the first UniProt was ever passed
        to get_related_by_uniprot_split / get_related_by_uniprot."""
        entry_data = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"diffrn_resolution_high": {"value": 2.0}, "nonpolymer_bound_components": []},
            "refine": [],
            "pdbx_vrpt_summary_geometry": [],
            "pdbx_vrpt_summary_diffraction": [],
            "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1", "2"], "non_polymer_entity_ids": []},
        }
        def _entity(uniprot_id, seq="A" * 100):
            return {
                "entity_poly": {"rcsb_entity_polymer_type": "Protein", "pdbx_seq_one_letter_code_can": seq},
                "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [uniprot_id], "bird_id": None},
                "rcsb_target_cofactors": [],
                "rcsb_polymer_entity_feature": [],
            }
        client = MagicMock()
        client.get.side_effect = [entry_data, _entity("P11111"), _entity("P22222"), {"1abc": []}]
        client.post.return_value = {"result_set": []}
        enrich_row(
            row={"PDBID": "1ABC", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        # Extract UniProt values searched across all POST calls
        searched_uniprots = set()
        for call in client.post.call_args_list:
            payload = call.args[1] if len(call.args) > 1 else call.kwargs.get("json", {})
            payload_str = str(payload)
            if "P11111" in payload_str:
                searched_uniprots.add("P11111")
            if "P22222" in payload_str:
                searched_uniprots.add("P22222")
        assert "P11111" in searched_uniprots, "first resolved UniProt must be used in related search"
        assert "P22222" in searched_uniprots, "second resolved UniProt must be used in related search — bug: was skipped"

    def test_uniprot_resolved_from_entity_when_missing(self):
        entry_data = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"diffrn_resolution_high": {"value": 2.0}, "nonpolymer_bound_components": []},
            "refine": [],
            "pdbx_vrpt_summary_geometry": [],
            "pdbx_vrpt_summary_diffraction": [],
            "rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"], "non_polymer_entity_ids": []},
        }
        entity_data = {
            "entity_poly": {"rcsb_entity_polymer_type": "Protein", "pdbx_seq_one_letter_code_can": "A" * 100},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["P99999"], "bird_id": None},
            "rcsb_target_cofactors": [],
            "rcsb_polymer_entity_feature": [],
        }
        client = MagicMock()
        client.get.side_effect = [entry_data, entity_data, {f"1abc": []}]
        client.post.return_value = {"result_set": []}
        result = enrich_row(
            row={"PDBID": "1ABC", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert "P99999" in result["Uniprot"]

    def test_ligand_detail_cols_are_none_on_primary(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        for col in LIGAND_DETAIL_COLS:
            assert result[col] is None, f"Expected {col} to be None on primary row"

    def test_related_tag_cols_are_none_on_primary(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert result["related_pdb_ids"] is None
        assert result["fulllength_pdb_ids"] is None

    def test_internal_related_ligand_entries_attached(self):
        client = self._mock_full_run()
        result = enrich_row(
            row=self._base_row(),
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=25,
        )
        assert "_sibling_ligand_entries" in result
        assert "_fulllength_ligand_entries" in result


# ---------------------------------------------------------------------------
# _fetch_related_ligand_data
# ---------------------------------------------------------------------------

class TestFetchRelatedLigandData:
    def _entry_data(self, np_entity_ids=None, poly_entity_ids=None):
        return {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {
                "diffrn_resolution_high": {"value": 2.0},
                "nonpolymer_bound_components": ["LIG"] if np_entity_ids else [],
            },
            "refine": [],
            "pdbx_vrpt_summary_geometry": [],
            "pdbx_vrpt_summary_diffraction": [],
            "rcsb_entry_container_identifiers": {
                "polymer_entity_ids": poly_entity_ids or [],
                "non_polymer_entity_ids": np_entity_ids or [],
            },
        }

    def test_404_returns_no_ligands(self):
        client = MagicMock()
        client.get.return_value = None
        result = _fetch_related_ligand_data(client, "XXXX")
        assert result["has_ligands"] is False
        assert result["ligand_metrics"] == []
        assert result["peptide_entities"] == []
        assert result["entry_quality"] == {}

    def test_entry_with_no_nonpolymer_entities(self):
        client = MagicMock()
        client.get.return_value = self._entry_data()
        result = _fetch_related_ligand_data(client, "1ABC")
        assert result["has_ligands"] is False

    def test_entry_quality_populated(self):
        client = MagicMock()
        client.get.return_value = self._entry_data()
        result = _fetch_related_ligand_data(client, "1ABC")
        assert "exp_method" in result["entry_quality"]
        assert "resolution_A" in result["entry_quality"]
        assert result["entry_quality"]["exp_method"] == "X-RAY DIFFRACTION"
        assert result["entry_quality"]["resolution_A"] == 2.0

    def test_entry_with_interesting_ligand_with_contacts(self):
        entry = self._entry_data(np_entity_ids=["1"])
        instance = {
            "rcsb_nonpolymer_entity_container_identifiers": {
                "nonpolymer_comp_id": "LIG", "asym_ids": ["A"]
            }
        }
        inst_data = {
            "rcsb_nonpolymer_instance_validation_score": [{
                "RSCC": 0.9, "RSR": 0.1, "mogul_bonds_RMSZ": 1.0,
                "mogul_angles_RMSZ": 1.0, "intermolecular_clashes": 0,
                "is_subject_of_investigation": True,
            }],
            "rcsb_target_neighbors": [
                {"distance": 3.5, "target_asym_id": "B", "target_seq_id": 10,
                 "target_comp_id": "ALA", "target_auth_seq_id": "10"}
            ],
        }
        chain_inst = {"rcsb_polymer_instance_feature": []}
        client = MagicMock()
        client.get.side_effect = [entry, instance, inst_data, chain_inst]
        result = _fetch_related_ligand_data(client, "1ABC")
        assert result["has_ligands"] is True
        assert len(result["ligand_metrics"]) == 1

    def test_entry_with_non_interesting_ligand_only(self):
        # A ligand marked is_interesting=False and contact_residue_count=0 → has_ligands=False
        entry = self._entry_data(np_entity_ids=["1"])
        instance = {
            "rcsb_nonpolymer_entity_container_identifiers": {
                "nonpolymer_comp_id": "HOH", "asym_ids": ["W"]
            }
        }
        inst_data = {
            "rcsb_nonpolymer_instance_validation_score": [{
                "RSCC": 0.9, "RSR": 0.1, "mogul_bonds_RMSZ": 1.0,
                "mogul_angles_RMSZ": 1.0, "intermolecular_clashes": 0,
                "is_subject_of_investigation": False,
            }],
            "rcsb_target_neighbors": [],
        }
        client = MagicMock()
        client.get.side_effect = [entry, instance, inst_data]
        result = _fetch_related_ligand_data(client, "1ABC")
        assert result["has_ligands"] is False

    def test_pdb_id_preserved(self):
        client = MagicMock()
        client.get.return_value = None
        result = _fetch_related_ligand_data(client, "4XYZ")
        assert result["pdb_id"] == "4XYZ"


# ---------------------------------------------------------------------------
# build_ligand_rows tags parameter
# ---------------------------------------------------------------------------

class TestBuildLigandRowsTags:
    def _all_cols(self):
        from rcsb_enrichment.cli import _AUGMENTED_COLS
        return ("pdb_id", "uniprot") + _AUGMENTED_COLS

    def _good_metric(self):
        return {
            "ligand_id": "LIG", "chain_id": "A", "is_interesting": True,
            "contact_residue_count": 3, "binding_quality": "good",
            "ligand_rscc": 0.90, "ligand_rsr": 0.12, "ligand_rmsz_bonds": 1.0,
            "ligand_rmsz_angles": 1.0, "ligand_intermolecular_clashes": 0,
            "contact_outlier_fraction": 0.05, "contact_residues": "ALA10(A)",
        }

    def test_tags_applied_to_small_molecule_row(self):
        rows = build_ligand_rows(
            "1ABC", [self._good_metric()], [], self._all_cols(),
            tags={"related_pdb_ids": "4SIB"}
        )
        assert rows[0]["related_pdb_ids"] == "4SIB"
        assert rows[0]["fulllength_pdb_ids"] is None

    def test_crystal_quality_tags_propagated(self):
        rows = build_ligand_rows(
            "1ABC", [self._good_metric()], [], self._all_cols(),
            tags={"related_pdb_ids": "4SIB", "exp_method": "X-RAY DIFFRACTION",
                  "resolution_A": 1.8, "r_work": 0.19, "r_free": 0.23,
                  "clashscore": 2.1, "rsrz_outliers_pct": 3.0}
        )
        assert rows[0]["exp_method"] == "X-RAY DIFFRACTION"
        assert rows[0]["resolution_A"] == pytest.approx(1.8)
        assert rows[0]["clashscore"] == pytest.approx(2.1)

    def test_tags_applied_to_peptide_row(self):
        rows = build_ligand_rows(
            "1ABC", [], [{"bird_id": "PRD_000001", "sequence": "ACDE"}], self._all_cols(),
            tags={"fulllength_pdb_ids": "8FUL"}
        )
        assert rows[0]["fulllength_pdb_ids"] == "8FUL"
        assert rows[0]["related_pdb_ids"] is None

    def test_no_tags_leaves_columns_none(self):
        rows = build_ligand_rows(
            "1ABC", [self._good_metric()], [], self._all_cols()
        )
        assert rows[0]["related_pdb_ids"] is None
        assert rows[0]["fulllength_pdb_ids"] is None

    def test_parent_pdb_id_is_input_pdb_not_related(self):
        rows = build_ligand_rows(
            "1ABC", [self._good_metric()], [], self._all_cols(),
            tags={"related_pdb_ids": "4SIB"}
        )
        assert rows[0]["parent_pdb_id"] == "1ABC"


# ---------------------------------------------------------------------------
# _entity_matches_names
# ---------------------------------------------------------------------------

class TestEntityMatchesNames:
    """Unit tests for the whitespace-token entity name filter.

    The filter is intentionally conservative: 'Tubulin' must not match
    'Tubulin-Tyrosine Ligase' because hyphenated words are single tokens.
    """

    def test_exact_word_matches(self):
        assert _entity_matches_names("Tubulin alpha-1B chain", ["Tubulin"])

    def test_case_insensitive(self):
        assert _entity_matches_names("Tubulin alpha-1B chain", ["tubulin"])
        assert _entity_matches_names("tubulin alpha-1b chain", ["Tubulin"])

    def test_hyphenated_compound_does_not_match_prefix(self):
        # 'Tubulin-Tyrosine' is one token — 'Tubulin' alone must NOT match
        assert not _entity_matches_names("Tubulin-Tyrosine Ligase", ["Tubulin"])

    def test_unrelated_description_does_not_match(self):
        assert not _entity_matches_names("Stathmin-4", ["Tubulin"])

    def test_multiple_filters_any_match_is_sufficient(self):
        assert _entity_matches_names("Stathmin-4", ["Tubulin", "Stathmin-4"])

    def test_empty_filters_list_returns_false(self):
        assert not _entity_matches_names("Tubulin alpha chain", [])

    def test_empty_description_returns_false(self):
        assert not _entity_matches_names("", ["Tubulin"])

    def test_5s5v_like_descriptions(self):
        """Replicate the four 5S5V entity descriptions; filter 'Tubulin' → exactly 2 matches."""
        descriptions = [
            "Tubulin alpha-1B chain",
            "Tubulin beta-2B chain",
            "Stathmin-4",
            "Tubulin-Tyrosine Ligase",
        ]
        matched = [d for d in descriptions if _entity_matches_names(d, ["Tubulin"])]
        assert matched == ["Tubulin alpha-1B chain", "Tubulin beta-2B chain"]


# ---------------------------------------------------------------------------
# enrich_row entity_name_filters integration
# ---------------------------------------------------------------------------

class TestEnrichRowEntityNameFilters:
    """5S5V has four protein entities (alpha-tubulin, beta-tubulin, Stathmin-4,
    Tubulin-Tyrosine Ligase).  With filter ['Tubulin'] only the two tubulin
    chains should be retained as receptor entities, affecting which UniProt IDs
    are collected and which are searched for related entries.

    All API calls are mocked; mock data mirrors the real 5S5V entity layout.
    """

    # CSV used in this test: data/examples/test_entity_filter.csv
    #   PDBID,Uniprot
    #   5S5V,

    def _entry_data(self):
        return {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"diffrn_resolution_high": {"value": 2.3}, "nonpolymer_bound_components": []},
            "refine": [{"ls_R_factor_R_work": 0.19, "ls_R_factor_R_free": 0.23}],
            "pdbx_vrpt_summary_geometry": [
                {"clashscore": 3.0, "percent_ramachandran_outliers": 0.2,
                 "percent_rotamer_outliers": 1.0, "bonds_RMSZ": 0.9, "angles_RMSZ": 1.1}
            ],
            "pdbx_vrpt_summary_diffraction": [{"percent_RSRZ_outliers": 4.0}],
            "rcsb_entry_container_identifiers": {
                "polymer_entity_ids": ["1", "2", "3", "4"],
                "non_polymer_entity_ids": [],
            },
        }

    def _entity(self, uniprot_id, description, seq="A" * 200):
        return {
            "entity_poly": {"rcsb_entity_polymer_type": "Protein", "pdbx_seq_one_letter_code_can": seq},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [uniprot_id], "bird_id": None},
            "rcsb_polymer_entity": {"pdbx_description": description},
            "rcsb_entity_source_organism": [],
            "rcsb_target_cofactors": [],
            "rcsb_polymer_entity_feature": [],
        }

    def _make_client(self):
        client = MagicMock()
        client.get.side_effect = [
            self._entry_data(),
            self._entity("P81947", "Tubulin alpha-1B chain"),   # entity 1 — kept
            self._entity("Q6B856", "Tubulin beta-2B chain"),    # entity 2 — kept
            self._entity("P63043", "Stathmin-4"),               # entity 3 — filtered out
            self._entity("E1BQ43", "Tubulin-Tyrosine Ligase"),  # entity 4 — filtered out
            {"5s5v": []},                                       # PDBe binding sites
        ]
        client.post.return_value = {"result_set": []}
        return client

    def test_filter_tubulin_retains_two_entities(self):
        client = self._make_client()
        result = enrich_row(
            row={"PDBID": "5S5V", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=5,
            entity_name_filters=["Tubulin"],
        )
        resolved = result["Uniprot"]
        assert "P81947" in resolved, "alpha-tubulin UniProt must be resolved"
        assert "Q6B856" in resolved, "beta-tubulin UniProt must be resolved"
        assert "P63043" not in resolved, "Stathmin-4 must be filtered out"
        assert "E1BQ43" not in resolved, "Tubulin-Tyrosine Ligase must be filtered out"

    def test_no_filter_retains_all_four_entities(self):
        client = self._make_client()
        result = enrich_row(
            row={"PDBID": "5S5V", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=5,
            entity_name_filters=None,
        )
        resolved = result["Uniprot"]
        assert "P81947" in resolved
        assert "Q6B856" in resolved
        assert "P63043" in resolved
        assert "E1BQ43" in resolved

    def test_filter_only_searches_matching_uniprots(self):
        """Only the two tubulin UniProt IDs should appear in any search API call."""
        client = self._make_client()
        enrich_row(
            row={"PDBID": "5S5V", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=5,
            entity_name_filters=["Tubulin"],
        )
        searched = set()
        for call in client.post.call_args_list:
            payload_str = str(call.args[1] if len(call.args) > 1 else call.kwargs.get("json", {}))
            for uid in ("P81947", "Q6B856", "P63043", "E1BQ43"):
                if uid in payload_str:
                    searched.add(uid)
        assert "P81947" in searched or "Q6B856" in searched, "tubulin IDs must be searched"
        assert "P63043" not in searched, "Stathmin-4 must not be searched"
        assert "E1BQ43" not in searched, "Tubulin-Tyrosine Ligase must not be searched"


# ---------------------------------------------------------------------------
# _collect_entity_names
# ---------------------------------------------------------------------------

class TestCollectEntityNames:
    """Unit tests for the entity-name aggregation helper.

    Covers polymer receptors, peptide ligands, non-polymer ligands, deduplication,
    and ordering (polymer first, then non-polymer).
    """

    def test_polymer_descriptions_included(self):
        entities = [{"description": "Tubulin alpha-1B chain"}, {"description": "Tubulin beta-2B chain"}]
        assert _collect_entity_names(entities, [], []) == "Tubulin alpha-1B chain,Tubulin beta-2B chain"

    def test_peptide_descriptions_included(self):
        peptides = [{"description": "Short peptide inhibitor"}]
        assert _collect_entity_names([], peptides, []) == "Short peptide inhibitor"

    def test_nonpolymer_descriptions_included(self):
        metrics = [{"description": "GUANOSINE-5'-TRIPHOSPHATE"}, {"description": "MAGNESIUM ION"}]
        assert _collect_entity_names([], [], metrics) == "GUANOSINE-5'-TRIPHOSPHATE,MAGNESIUM ION"

    def test_polymer_before_nonpolymer(self):
        entities = [{"description": "HIV-1 Protease"}]
        metrics = [{"description": "Indinavir"}]
        result = _collect_entity_names(entities, [], metrics)
        assert result.index("HIV-1 Protease") < result.index("Indinavir")

    def test_duplicates_deduplicated(self):
        # same description across polymer chains (e.g. homodimer)
        entities = [{"description": "Kinase domain"}, {"description": "Kinase domain"}]
        assert _collect_entity_names(entities, [], []) == "Kinase domain"

    def test_empty_descriptions_skipped(self):
        entities = [{"description": ""}, {"description": "Kinase domain"}]
        assert _collect_entity_names(entities, [], []) == "Kinase domain"

    def test_missing_description_key_skipped(self):
        entities = [{"uniprot_ids": ["P12345"]}, {"description": "Kinase domain"}]
        assert _collect_entity_names(entities, [], []) == "Kinase domain"

    def test_all_empty_returns_empty_string(self):
        assert _collect_entity_names([], [], []) == ""

    def test_5s5v_like_four_entities(self):
        """All four 5S5V polymer entities; no non-polymer."""
        entities = [
            {"description": "Tubulin alpha-1B chain"},
            {"description": "Tubulin beta-2B chain"},
            {"description": "Stathmin-4"},
            {"description": "Tubulin-Tyrosine Ligase"},
        ]
        result = _collect_entity_names(entities, [], [])
        assert result == (
            "Tubulin alpha-1B chain,Tubulin beta-2B chain,"
            "Stathmin-4,Tubulin-Tyrosine Ligase"
        )


# ---------------------------------------------------------------------------
# entity_names propagation through enrich_row and related-entry sub-rows
# ---------------------------------------------------------------------------

class TestEntityNamesInEnrichRow:
    """Verify entity_names is set on primary rows from enrich_row, and that
    _fetch_related_ligand_data returns it so cli.py can propagate it to
    sibling/full-length sub-rows.

    Test CSV used: data/examples/test_entity_filter.csv
        PDBID,Uniprot
        5S5V,
    """

    def _entry_data(self, polymer_ids, np_ids=None):
        return {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"diffrn_resolution_high": {"value": 2.0}, "nonpolymer_bound_components": []},
            "refine": [],
            "pdbx_vrpt_summary_geometry": [],
            "pdbx_vrpt_summary_diffraction": [],
            "rcsb_entry_container_identifiers": {
                "polymer_entity_ids": polymer_ids,
                "non_polymer_entity_ids": np_ids or [],
            },
        }

    def _poly_entity(self, uniprot_id, description, seq="A" * 100):
        return {
            "entity_poly": {"rcsb_entity_polymer_type": "Protein", "pdbx_seq_one_letter_code_can": seq},
            "rcsb_polymer_entity_container_identifiers": {"uniprot_ids": [uniprot_id], "bird_id": None},
            "rcsb_polymer_entity": {"pdbx_description": description},
            "rcsb_entity_source_organism": [],
            "rcsb_target_cofactors": [],
            "rcsb_polymer_entity_feature": [],
        }

    def test_entity_names_on_primary_row(self):
        client = MagicMock()
        client.get.side_effect = [
            self._entry_data(["1", "2"]),
            self._poly_entity("P81947", "Tubulin alpha-1B chain"),
            self._poly_entity("Q6B856", "Tubulin beta-2B chain"),
            {"5s5v": []},  # PDBe
        ]
        client.post.return_value = {"result_set": []}
        result = enrich_row(
            row={"PDBID": "5S5V", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=5,
        )
        assert result["entity_names"] == "Tubulin alpha-1B chain,Tubulin beta-2B chain"

    def test_entity_names_empty_on_404(self):
        client = MagicMock()
        client.get.return_value = None
        client.post.return_value = {"result_set": []}
        result = enrich_row(
            row={"PDBID": "1ABC", "Uniprot": ""},
            client=client,
            pdb_col="PDBID",
            uniprot_col="Uniprot",
            seq_identity=0.9,
            max_related=5,
        )
        assert result["entity_names"] == ""

    def test_fetch_related_ligand_data_returns_entity_names(self):
        """_fetch_related_ligand_data must include entity_names in its return dict."""
        client = MagicMock()
        client.get.side_effect = [
            self._entry_data(["1"], np_ids=[]),
            self._poly_entity("P12345", "HIV-1 Protease"),
        ]
        data = _fetch_related_ligand_data(client, "1HSG")
        assert "entity_names" in data
        assert "HIV-1 Protease" in data["entity_names"]
