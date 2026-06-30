"""Tests for entities.get_polymer_entities and extract_direct_binders."""

import pytest
from unittest.mock import MagicMock
from rcsb_enrichment.entities import _PEPTIDE_LEN_THRESHOLD, extract_direct_binders, get_polymer_entities


def _poly_entity(
    uniprot_ids=None,
    sequence="ACDEFGHIKLMNPQRSTVWY",
    bird_id=None,
    cofactors=None,
    polymer_type="Protein",
):
    data = {
        "entity_poly": {
            "rcsb_entity_polymer_type": polymer_type,
            "pdbx_seq_one_letter_code_can": sequence,
        },
        "rcsb_polymer_entity_container_identifiers": {
            "uniprot_ids": uniprot_ids or [],
            "bird_id": bird_id,
        },
        "rcsb_target_cofactors": cofactors or [],
        "rcsb_polymer_entity_feature": [],
    }
    return data


class TestGetPolymerEntities:
    def test_normal_protein_entity(self):
        client = MagicMock()
        client.get.return_value = _poly_entity(uniprot_ids=["P12345"], sequence="A" * 100)
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert len(result) == 1
        assert result[0]["uniprot_ids"] == ["P12345"]
        assert result[0]["is_peptide_ligand"] is False

    def test_non_protein_skipped(self):
        client = MagicMock()
        client.get.return_value = _poly_entity(polymer_type="RNA")
        assert get_polymer_entities(client, "1ABC", ["1"]) == []

    def test_bird_annotated_is_peptide(self):
        client = MagicMock()
        client.get.return_value = _poly_entity(bird_id="PRD_000001", sequence="ACDEF")
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert result[0]["is_peptide_ligand"] is True
        assert result[0]["bird_id"] == "PRD_000001"

    def test_short_no_uniprot_is_peptide(self):
        client = MagicMock()
        short_seq = "A" * _PEPTIDE_LEN_THRESHOLD
        client.get.return_value = _poly_entity(uniprot_ids=[], sequence=short_seq)
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert result[0]["is_peptide_ligand"] is True

    def test_short_with_uniprot_is_not_peptide(self):
        client = MagicMock()
        short_seq = "A" * _PEPTIDE_LEN_THRESHOLD
        client.get.return_value = _poly_entity(uniprot_ids=["P99999"], sequence=short_seq)
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert result[0]["is_peptide_ligand"] is False

    def test_just_over_threshold_is_not_peptide(self):
        client = MagicMock()
        seq = "A" * (_PEPTIDE_LEN_THRESHOLD + 1)
        client.get.return_value = _poly_entity(uniprot_ids=[], sequence=seq)
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert result[0]["is_peptide_ligand"] is False

    def test_cofactors_filtered_by_neighbor_flag(self):
        cofactors = [
            {"neighbor_flag": "N", "cofactor_name": "DirectBinder", "cofactor_SMILES": "C",
             "cofactor_InChIKey": "IK1", "cofactor_chem_comp_id": "LIG", "resource_name": "ChEMBL"},
            {"neighbor_flag": "Y", "cofactor_name": "Neighbour", "cofactor_SMILES": "CC",
             "cofactor_InChIKey": "IK2", "cofactor_chem_comp_id": "NEI", "resource_name": "DrugBank"},
        ]
        client = MagicMock()
        client.get.return_value = _poly_entity(uniprot_ids=["P12345"], sequence="A" * 200, cofactors=cofactors)
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert len(result[0]["cofactors"]) == 1
        assert result[0]["cofactors"][0]["name"] == "DirectBinder"

    def test_site_features_extracted(self):
        data = _poly_entity(uniprot_ids=["P12345"], sequence="A" * 200)
        data["rcsb_polymer_entity_feature"] = [
            {"type": "BINDING_SITE", "name": "ATP site", "feature_id": "feat1", "provenance_source": "UniProt"},
            {"type": "OTHER", "name": "ignored", "feature_id": "feat2", "provenance_source": "X"},
        ]
        client = MagicMock()
        client.get.return_value = data
        result = get_polymer_entities(client, "1ABC", ["1"])
        assert len(result[0]["site_features"]) == 1
        assert result[0]["site_features"][0]["name"] == "ATP site"

    def test_peptide_entity_skips_cofactor_fetch(self):
        client = MagicMock()
        client.get.return_value = _poly_entity(bird_id="PRD_000001", sequence="ACDEF")
        result = get_polymer_entities(client, "1ABC", ["1"])
        # Peptide entities don't have cofactors or site_features keys
        assert "cofactors" not in result[0]
        assert "site_features" not in result[0]

    def test_404_entity_skipped(self):
        client = MagicMock()
        client.get.return_value = None
        assert get_polymer_entities(client, "1ABC", ["1"]) == []

    def test_multiple_entities(self):
        client = MagicMock()
        client.get.side_effect = [
            _poly_entity(uniprot_ids=["P1"], sequence="A" * 200),
            _poly_entity(uniprot_ids=["P2"], sequence="A" * 150),
        ]
        result = get_polymer_entities(client, "1ABC", ["1", "2"])
        assert len(result) == 2


class TestExtractDirectBinders:
    def test_returns_direct_binders(self):
        entry_data = {"rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"]}}
        entity_data = _poly_entity(
            uniprot_ids=["P12345"],
            sequence="A" * 200,
            cofactors=[{
                "neighbor_flag": "N",
                "cofactor_name": "Drug",
                "cofactor_SMILES": "C",
                "cofactor_InChIKey": "IK1",
                "chem_comp_id": "DRG",
                "resource_name": "ChEMBL",
            }],
        )
        client = MagicMock()
        client.get.side_effect = [entry_data, entity_data]
        binders = extract_direct_binders(client, "1ABC")
        assert len(binders) == 1
        assert binders[0]["name"] == "Drug"
        assert binders[0]["from_pdb"] == "1ABC"

    def test_neighbour_flag_y_excluded(self):
        entry_data = {"rcsb_entry_container_identifiers": {"polymer_entity_ids": ["1"]}}
        entity_data = _poly_entity(
            uniprot_ids=["P12345"],
            sequence="A" * 200,
            cofactors=[{"neighbor_flag": "Y", "cofactor_name": "Neighbour",
                        "cofactor_SMILES": "C", "cofactor_InChIKey": "IK",
                        "chem_comp_id": "NEI", "resource_name": "DrugBank"}],
        )
        client = MagicMock()
        client.get.side_effect = [entry_data, entity_data]
        assert extract_direct_binders(client, "1ABC") == []

    def test_404_entry_returns_empty(self):
        client = MagicMock()
        client.get.return_value = None
        assert extract_direct_binders(client, "1ABC") == []
