"""Tests for cli.detect_pdb_col, detect_uniprot_col, and main."""

import csv
import io
import pytest
from unittest.mock import MagicMock, patch
from rcsb_enrichment.cli import (
    detect_pdb_col,
    detect_uniprot_col,
)


# ---------------------------------------------------------------------------
# detect_pdb_col
# ---------------------------------------------------------------------------

class TestDetectPdbCol:
    def test_hint_found(self):
        assert detect_pdb_col(["pdb_id", "name"], "pdb_id") == "pdb_id"

    def test_hint_not_found_raises(self):
        with pytest.raises(ValueError, match="pdb_id"):
            detect_pdb_col(["name", "value"], "pdb_id")

    def test_auto_detect_pdb(self):
        assert detect_pdb_col(["name", "pdb", "value"], None) == "pdb"

    def test_auto_detect_pdbid(self):
        assert detect_pdb_col(["PDBID", "name"], None) == "PDBID"

    def test_auto_detect_case_insensitive(self):
        assert detect_pdb_col(["PDB_ID", "name"], None) == "PDB_ID"

    def test_auto_detect_fails_raises(self):
        with pytest.raises(ValueError, match="--pdb-col"):
            detect_pdb_col(["protein", "mw"], None)


# ---------------------------------------------------------------------------
# detect_uniprot_col
# ---------------------------------------------------------------------------

class TestDetectUniprotCol:
    def test_hint_found(self):
        assert detect_uniprot_col(["uniprot", "name"], "uniprot") == "uniprot"

    def test_hint_not_found_returns_none(self):
        assert detect_uniprot_col(["name"], "uniprot") is None

    def test_auto_detect_uniprot_id(self):
        assert detect_uniprot_col(["protein", "uniprot_id"], None) == "uniprot_id"

    def test_auto_detect_accession(self):
        assert detect_uniprot_col(["accession", "name"], None) == "accession"

    def test_no_match_returns_none(self):
        assert detect_uniprot_col(["protein", "pdb"], None) is None



# ---------------------------------------------------------------------------
# main — end-to-end with mocked API
# ---------------------------------------------------------------------------

class TestMain:
    def _minimal_enriched_row(self, pdb_id, uniprot):
        """Returns what enrich_row produces for a minimal valid entry."""
        from rcsb_enrichment.cli import _AUGMENTED_COLS
        row = {"PDBID": pdb_id, "Uniprot": uniprot}
        for col in _AUGMENTED_COLS:
            row[col] = None
        row["row_type"] = "primary"
        row["exp_method"] = "X-RAY DIFFRACTION"
        row["resolution_A"] = 2.0
        row["_seq_len"] = 100
        row["_resolved_uniprot_ids"] = [uniprot]
        row["_ligand_metrics"] = []
        row["_peptide_entities"] = []
        row["_sibling_ligand_entries"] = []
        row["_fulllength_ligand_entries"] = []
        return row

    def test_main_writes_csv(self, tmp_path):
        import sys
        input_csv = tmp_path / "input.csv"
        output_csv = tmp_path / "output.csv"
        input_csv.write_text("PDBID,Uniprot\n1ABC,P12345\n")

        with patch("rcsb_enrichment.cli.RCSBClient") as MockClient, \
             patch("rcsb_enrichment.cli.enrich_row") as mock_enrich:

            mock_enrich.return_value = self._minimal_enriched_row("1ABC", "P12345")

            sys.argv = ["rcsb-enrich", "-i", str(input_csv), "-o", str(output_csv),
                        "--pdb-col", "PDBID", "--uniprot-col", "Uniprot"]
            from rcsb_enrichment.cli import main
            main()

        assert output_csv.exists()
        with open(output_csv) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) >= 1
        assert any(r.get("PDBID") == "1ABC" or r.get("parent_pdb_id") == "1ABC" for r in rows)

    def test_main_handles_empty_pdb_id(self, tmp_path):
        import sys
        input_csv = tmp_path / "input.csv"
        output_csv = tmp_path / "output.csv"
        # Row 1 has a blank PDBID, row 2 is valid
        input_csv.write_text("PDBID,Uniprot\n,P12345\n1ABC,P12345\n")

        with patch("rcsb_enrichment.cli.RCSBClient"), \
             patch("rcsb_enrichment.cli.enrich_row") as mock_enrich:

            mock_enrich.return_value = self._minimal_enriched_row("1ABC", "P12345")

            sys.argv = ["rcsb-enrich", "-i", str(input_csv), "-o", str(output_csv),
                        "--pdb-col", "PDBID", "--uniprot-col", "Uniprot"]
            from rcsb_enrichment.cli import main
            main()

        # pandas reads a blank PDBID cell as NaN → normalises to "NAN" → both rows reach enrich_row
        assert mock_enrich.call_count == 2
        assert output_csv.exists()
        with open(output_csv) as f:
            rows = list(csv.DictReader(f))
        assert len(rows) >= 1
