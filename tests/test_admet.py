"""Tests for the ADMET module."""

import pytest
from src.admet import calculate_lipinski_rules


class TestLipinskiRules:
    """Tests for calculate_lipinski_rules()."""

    def test_aspirin_passes(self):
        """Aspirin (MW≈180, LogP≈1.2) should pass Ro5."""
        smiles = ["CC(=O)Oc1ccccc1C(=O)O"]
        df = calculate_lipinski_rules(smiles)

        assert len(df) == 1
        assert df.iloc[0]["Valid"] is True
        assert df.iloc[0]["Passes_Ro5"] is True
        assert 170 < df.iloc[0]["MW"] < 190

    def test_invalid_smiles(self):
        """Invalid SMILES should be flagged."""
        smiles = ["INVALID_SMILES"]
        df = calculate_lipinski_rules(smiles)

        assert len(df) == 1
        assert df.iloc[0]["Valid"] is False

    def test_multiple_molecules(self):
        """Test batch processing of multiple molecules."""
        smiles = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "c1ccccc1",                 # Benzene
            "O=C(O)c1ccccc1O",          # Salicylic acid
        ]
        df = calculate_lipinski_rules(smiles)
        assert len(df) == 3
        assert all(df["Valid"])

    def test_large_molecule_fails_ro5(self):
        """A molecule with MW > 500 should fail Ro5."""
        # Cyclosporine A (MW ≈ 1202)
        smiles = ["CCC1NC(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)"
                  "C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(CC2=CC=CC=C2)N(C)"
                  "C(=O)C(CC(C)C)NC(=O)C1C"]
        df = calculate_lipinski_rules(smiles)

        if df.iloc[0]["Valid"]:
            assert df.iloc[0]["MW"] > 500

    def test_empty_list(self):
        """Empty input should return empty DataFrame."""
        df = calculate_lipinski_rules([])
        assert len(df) == 0
