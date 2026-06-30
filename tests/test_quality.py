"""Tests for quality.traffic_light, _parse_residue_features, get_entry_quality, get_ligand_quality."""

import pytest
from unittest.mock import MagicMock, call
from rcsb_enrichment.quality import (
    CONTACT_DIST_CUTOFF,
    _DENSE_FEATURE_TYPES,
    _SPARSE_FEATURE_TYPES,
    _get_contact_residue_outlier_fraction,
    _parse_residue_features,
    get_entry_quality,
    get_ligand_quality,
    traffic_light,
)


# ---------------------------------------------------------------------------
# traffic_light
# ---------------------------------------------------------------------------

class TestTrafficLight:
    def test_all_good(self):
        assert traffic_light(0.90, 0.15, 1.0, 1.0, 0, 0.05) == "good"

    def test_all_bad(self):
        assert traffic_light(0.40, 0.40, 3.0, 3.0, 5, 0.50) == "bad"

    def test_hard_override_rscc_below_0_50(self):
        # Even with otherwise good metrics, RSCC < 0.50 forces "bad"
        assert traffic_light(0.49, 0.10, 0.5, 0.5, 0, 0.0) == "bad"

    def test_rscc_exactly_0_50_not_overridden(self):
        # 0.50 is NOT below the threshold — should not hard-override
        result = traffic_light(0.50, 0.10, 0.5, 0.5, 0, 0.0)
        assert result in ("good", "fair")

    def test_fair_boundary(self):
        # Craft a case that lands between good and bad
        result = traffic_light(0.70, 0.28, 1.8, 1.8, 1, 0.15)
        assert result == "fair"

    def test_missing_rscc_excluded_from_mean(self):
        # Without RSCC the score is computed only from the remaining metrics
        result_with = traffic_light(0.95, 0.10, 1.0, 1.0, 0, 0.0)
        result_without = traffic_light(None, 0.10, 1.0, 1.0, 0, 0.0)
        assert result_with == result_without == "good"

    def test_all_none_except_contact_fraction(self):
        # contact_outlier_fraction always contributes; others None
        result = traffic_light(None, None, None, None, None, 0.05)
        assert result == "good"

    def test_clashes_boundary(self):
        # clashes=0 → score 0 (good); clashes=2 → score 1 (fair); clashes=3 → score 2 (bad)
        # With only one metric contributing the mean equals that metric's score:
        # mean 0 < 0.67 → good; mean 1 → fair; mean 2 ≥ 1.33 → bad
        assert traffic_light(None, None, None, None, 0, 0.0) == "good"
        # contact_outlier_fraction=0.0 also contributes score 0 → mean = (0+1)/2 = 0.5 → good
        assert traffic_light(None, None, None, None, 2, 0.0) == "good"
        # clashes=3 → score 2; contact_frac=0 → score 0; mean = 1.0 → fair
        assert traffic_light(None, None, None, None, 3, 0.0) == "fair"

    def test_contact_fraction_thresholds(self):
        assert traffic_light(None, None, None, None, None, 0.10) == "good"
        assert traffic_light(None, None, None, None, None, 0.25) in ("good", "fair")
        assert traffic_light(None, None, None, None, None, 0.30) == "bad"


# ---------------------------------------------------------------------------
# _parse_residue_features
# ---------------------------------------------------------------------------

class TestParseResidueFeatures:
    def test_dense_feature(self):
        features = [{
            "type": "RSRZ",
            "feature_positions": [{"beg_seq_id": 10, "values": [0.5, 1.2, 0.3]}],
        }]
        result = _parse_residue_features(features)
        assert result[10]["RSRZ"] == 0.5
        assert result[11]["RSRZ"] == 1.2
        assert result[12]["RSRZ"] == 0.3

    def test_sparse_feature(self):
        features = [{
            "type": "RAMACHANDRAN_OUTLIER",
            "feature_positions": [{"beg_seq_id": 42, "values": [1]}],
        }]
        result = _parse_residue_features(features)
        assert result[42]["RAMACHANDRAN_OUTLIER"] == 1

    def test_multiple_feature_types_merged(self):
        features = [
            {"type": "RSRZ", "feature_positions": [{"beg_seq_id": 5, "values": [0.8]}]},
            {"type": "ROTAMER_OUTLIER", "feature_positions": [{"beg_seq_id": 5, "values": [1]}]},
        ]
        result = _parse_residue_features(features)
        assert "RSRZ" in result[5]
        assert "ROTAMER_OUTLIER" in result[5]

    def test_empty_values_skipped(self):
        features = [{"type": "RSRZ", "feature_positions": [{"beg_seq_id": 1, "values": []}]}]
        assert _parse_residue_features(features) == {}

    def test_missing_beg_seq_id_skipped(self):
        features = [{"type": "RSRZ", "feature_positions": [{"values": [0.5]}]}]
        assert _parse_residue_features(features) == {}

    def test_unknown_feature_type_ignored(self):
        features = [{"type": "UNKNOWN_TYPE", "feature_positions": [{"beg_seq_id": 1, "values": [1.0]}]}]
        assert _parse_residue_features(features) == {}

    def test_empty_input(self):
        assert _parse_residue_features([]) == {}


# ---------------------------------------------------------------------------
# _get_contact_residue_outlier_fraction
# ---------------------------------------------------------------------------

class TestContactResidueOutlierFraction:
    def _make_client(self, instance_data):
        client = MagicMock()
        client.get.return_value = instance_data
        return client

    def test_no_neighbors_returns_zero(self):
        client = self._make_client({})
        assert _get_contact_residue_outlier_fraction(client, "1ABC", [], {}) == 0.0

    def test_neighbors_beyond_cutoff_ignored(self):
        neighbors = [{"distance": CONTACT_DIST_CUTOFF + 0.1, "target_asym_id": "A", "target_seq_id": 10}]
        client = self._make_client({"rcsb_polymer_instance_feature": []})
        result = _get_contact_residue_outlier_fraction(client, "1ABC", neighbors, {})
        assert result == 0.0

    def test_contact_residue_with_outlier(self):
        neighbors = [{"distance": 3.5, "target_asym_id": "A", "target_seq_id": 10}]
        instance_data = {
            "rcsb_polymer_instance_feature": [{
                "type": "RAMACHANDRAN_OUTLIER",
                "feature_positions": [{"beg_seq_id": 10, "values": [1]}],
            }]
        }
        client = self._make_client(instance_data)
        result = _get_contact_residue_outlier_fraction(client, "1ABC", neighbors, {})
        assert result == 1.0

    def test_contact_residue_without_outlier(self):
        neighbors = [{"distance": 3.5, "target_asym_id": "A", "target_seq_id": 10}]
        instance_data = {"rcsb_polymer_instance_feature": []}
        client = self._make_client(instance_data)
        result = _get_contact_residue_outlier_fraction(client, "1ABC", neighbors, {})
        assert result == 0.0

    def test_cache_prevents_refetch(self):
        neighbors = [
            {"distance": 3.5, "target_asym_id": "A", "target_seq_id": 10},
            {"distance": 3.5, "target_asym_id": "A", "target_seq_id": 11},
        ]
        cache = {"A": {10: {"RAMACHANDRAN_OUTLIER": 1}, 11: {}}}
        client = MagicMock()
        result = _get_contact_residue_outlier_fraction(client, "1ABC", neighbors, cache)
        client.get.assert_not_called()
        assert result == pytest.approx(0.5)

    def test_fraction_calculation(self):
        neighbors = [
            {"distance": 3.0, "target_asym_id": "A", "target_seq_id": i}
            for i in range(1, 5)
        ]
        # residues 1 and 3 have outliers, 2 and 4 do not
        cache = {
            "A": {
                1: {"CLASHES": 1},
                2: {},
                3: {"ROTAMER_OUTLIER": 1},
                4: {},
            }
        }
        client = MagicMock()
        result = _get_contact_residue_outlier_fraction(client, "1ABC", neighbors, cache)
        assert result == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# get_entry_quality
# ---------------------------------------------------------------------------

class TestGetEntryQuality:
    def _entry_data(self, **overrides):
        base = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {
                "diffrn_resolution_high": {"value": 2.0},
                "nonpolymer_bound_components": ["ATP", "MG"],
            },
            "refine": [{"ls_R_factor_R_work": 0.18, "ls_R_factor_R_free": 0.22}],
            "pdbx_vrpt_summary_geometry": [{
                "clashscore": 2.1,
                "percent_ramachandran_outliers": 0.5,
                "percent_rotamer_outliers": 1.2,
                "bonds_RMSZ": 0.8,
                "angles_RMSZ": 1.1,
            }],
            "pdbx_vrpt_summary_diffraction": [{"percent_RSRZ_outliers": 3.5}],
            "rcsb_entry_container_identifiers": {
                "polymer_entity_ids": ["1", "2"],
                "non_polymer_entity_ids": ["3"],
            },
        }
        base.update(overrides)
        return base

    def test_full_entry(self):
        client = MagicMock()
        client.get.return_value = self._entry_data()
        result = get_entry_quality(client, "1ABC")
        assert result["exp_method"] == "X-RAY DIFFRACTION"
        assert result["resolution_A"] == 2.0
        assert result["r_work"] == pytest.approx(0.18)
        assert result["r_free"] == pytest.approx(0.22)
        assert result["clashscore"] == pytest.approx(2.1)
        assert result["rsrz_outliers_pct"] == pytest.approx(3.5)
        assert result["ligands_present"] == "ATP,MG"
        assert result["_polymer_entity_ids"] == ["1", "2"]
        assert result["_nonpolymer_entity_ids"] == ["3"]

    def test_returns_empty_dict_on_404(self):
        client = MagicMock()
        client.get.return_value = None
        assert get_entry_quality(client, "XXXX") == {}

    def test_missing_geometry_block_gives_none_fields(self):
        data = self._entry_data()
        data["pdbx_vrpt_summary_geometry"] = []
        client = MagicMock()
        client.get.return_value = data
        result = get_entry_quality(client, "1ABC")
        assert result["clashscore"] is None
        assert result["ramachandran_outliers_pct"] is None

    def test_rsrz_absent_for_non_xray(self):
        data = self._entry_data()
        data["pdbx_vrpt_summary_diffraction"] = []
        client = MagicMock()
        client.get.return_value = data
        result = get_entry_quality(client, "1ABC")
        assert result["rsrz_outliers_pct"] is None

    def test_no_ligands_gives_empty_string(self):
        data = self._entry_data()
        data["rcsb_entry_info"]["nonpolymer_bound_components"] = []
        client = MagicMock()
        client.get.return_value = data
        assert get_entry_quality(client, "1ABC")["ligands_present"] == ""

    def test_no_refine_block(self):
        data = self._entry_data()
        data["refine"] = []
        client = MagicMock()
        client.get.return_value = data
        result = get_entry_quality(client, "1ABC")
        assert result["r_work"] is None
        assert result["r_free"] is None


# ---------------------------------------------------------------------------
# get_ligand_quality
# ---------------------------------------------------------------------------

class TestGetLigandQuality:
    def _make_entity_data(self, comp_id="LIG", asym_ids=("A",)):
        return {
            "rcsb_nonpolymer_entity_container_identifiers": {
                "nonpolymer_comp_id": comp_id,
                "asym_ids": list(asym_ids),
            }
        }

    def _make_instance_data(self, rscc=0.88, rsr=0.12, isi=True, neighbors=None):
        return {
            "rcsb_nonpolymer_instance_validation_score": [{
                "RSCC": rscc,
                "RSR": rsr,
                "mogul_bonds_RMSZ": 1.0,
                "mogul_angles_RMSZ": 1.0,
                "intermolecular_clashes": 0,
                "is_subject_of_investigation": isi,
            }],
            "rcsb_target_neighbors": neighbors or [],
        }

    def test_basic_interesting_ligand(self):
        client = MagicMock()
        client.get.side_effect = [
            self._make_entity_data("LIG", ["A"]),
            self._make_instance_data(rscc=0.88, rsr=0.12, isi=True),
        ]
        results = get_ligand_quality(client, "1ABC", ["1"])
        assert len(results) == 1
        assert results[0]["ligand_id"] == "LIG"
        assert results[0]["is_interesting"] is True
        assert results[0]["ligand_rscc"] == pytest.approx(0.88)
        assert results[0]["binding_quality"] == "good"

    def test_non_interesting_flag(self):
        client = MagicMock()
        client.get.side_effect = [
            self._make_entity_data("HOH", ["W"]),
            self._make_instance_data(isi=False),
        ]
        results = get_ligand_quality(client, "1ABC", ["1"])
        assert results[0]["is_interesting"] is False

    def test_multiple_instances(self):
        client = MagicMock()
        client.get.side_effect = [
            self._make_entity_data("LIG", ["A", "B"]),
            self._make_instance_data(),
            self._make_instance_data(rscc=0.45),  # hard override → bad
        ]
        results = get_ligand_quality(client, "1ABC", ["1"])
        assert len(results) == 2
        assert results[0]["binding_quality"] == "good"
        assert results[1]["binding_quality"] == "bad"

    def test_404_entity_skipped(self):
        client = MagicMock()
        client.get.return_value = None
        results = get_ligand_quality(client, "1ABC", ["1"])
        assert results == []

    def test_contact_residues_formatted(self):
        neighbors = [{
            "distance": 3.5,
            "target_asym_id": "B",
            "target_seq_id": 10,
            "target_comp_id": "ALA",
            "target_auth_seq_id": "10",
        }]
        client = MagicMock()
        # entity fetch, then instance fetch, then chain feature fetch
        client.get.side_effect = [
            self._make_entity_data("LIG", ["A"]),
            self._make_instance_data(neighbors=neighbors),
            {"rcsb_polymer_instance_feature": []},
        ]
        results = get_ligand_quality(client, "1ABC", ["1"])
        assert results[0]["contact_residue_count"] == 1
        assert "ALA10(B)" in results[0]["contact_residues"]
