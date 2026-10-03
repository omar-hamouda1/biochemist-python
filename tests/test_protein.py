"""Tests for the protein module."""

import pytest
from pathlib import Path
from unittest.mock import patch

# We test import paths but actual PDB parsing needs Biopython + real files
from src.protein import load_protein, clean_protein, get_residue_count


PDB_FILE = Path("pdb/2zq2.pdb")


@pytest.mark.skipif(not PDB_FILE.exists(), reason="PDB file not available")
class TestLoadProtein:
    """Tests for load_protein()."""

    def test_load_returns_structure(self):
        structure = load_protein(PDB_FILE)
        assert structure is not None
        assert structure.id == "2zq2"

    def test_load_nonexistent_raises(self):
        with pytest.raises(FileNotFoundError):
            load_protein(Path("nonexistent.pdb"))

    def test_structure_has_models(self):
        structure = load_protein(PDB_FILE)
        models = list(structure.get_models())
        assert len(models) >= 1


@pytest.mark.skipif(not PDB_FILE.exists(), reason="PDB file not available")
class TestCleanProtein:
    """Tests for clean_protein()."""

    def test_clean_removes_water(self):
        structure = load_protein(PDB_FILE)
        cleaned = clean_protein(structure, keep_het=False)

        for model in cleaned:
            for chain in model:
                for residue in chain:
                    assert residue.id[0] != "W", "Water should be removed"

    def test_clean_keep_het(self):
        structure = load_protein(PDB_FILE)
        cleaned = clean_protein(structure, keep_het=True)
        # Should still have some residues
        assert get_residue_count(cleaned) > 0


@pytest.mark.skipif(not PDB_FILE.exists(), reason="PDB file not available")
class TestResidueCount:
    """Tests for get_residue_count()."""

    def test_count_positive(self):
        structure = load_protein(PDB_FILE)
        count = get_residue_count(structure)
        # Trypsin has ~223 residues
        assert count > 200
        assert count < 300
