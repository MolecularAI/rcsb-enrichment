"""Tests for ligand_filter.is_interesting_ligand."""

import pytest
from rcsb_enrichment.ligand_filter import _NON_INTERESTING_CCD, is_interesting_ligand


class TestIsInterestingLigand:
    def test_isi_true_overrides_exclusion_list(self):
        # Even a "non-interesting" CCD is interesting when ISI=True
        assert is_interesting_ligand("HOH", True) is True

    def test_isi_false_overrides_novel_ccd(self):
        assert is_interesting_ligand("NOVEL_DRUG", False) is False

    def test_isi_none_excluded_ccd(self):
        for ccd in ("HOH", "MG", "ATP", "GOL", "HEPES", "SO4"):
            assert is_interesting_ligand(ccd, None) is False, f"{ccd} should not be interesting"

    def test_isi_none_novel_ccd_is_interesting(self):
        assert is_interesting_ligand("XYZ", None) is True
        assert is_interesting_ligand("LIG", None) is True

    def test_case_insensitive_fallback(self):
        # Exclusion list check is case-insensitive via .upper()
        assert is_interesting_ligand("hoh", None) is False
        assert is_interesting_ligand("atp", None) is False

    def test_non_interesting_set_is_not_empty(self):
        assert len(_NON_INTERESTING_CCD) > 50

    def test_drug_like_ccd_not_in_exclusion_list(self):
        # Common drug-like CCDs that should never appear in the exclusion set
        for ccd in ("STI", "IMA", "LIG", "INH"):
            assert ccd not in _NON_INTERESTING_CCD
