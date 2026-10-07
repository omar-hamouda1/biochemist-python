"""Tests for the docking module."""

from pathlib import Path

import pytest

from src.docking import parse_smina_output, get_top_hits
from scripts.run_standardized_docking import parse_affinity


AFFINITIES_CSV = Path("docking/results/standardized_affinities.csv")

class TestParseSminaOutput:
    """Tests for parse_smina_output()."""

    def test_empty_file_returns_empty_list(self, tmp_path):
        empty_file = tmp_path / "empty.pdbqt"
        empty_file.write_text("")
        assert parse_smina_output(str(empty_file)) == []

    def test_nonexistent_file_returns_empty_list(self):
        assert parse_smina_output("nonexistent.pdbqt") == []

    def test_parses_valid_pdbqt(self, tmp_path):
        pdbqt = tmp_path / "test_out.pdbqt"
        pdbqt.write_text(
            "REMARK minimizedAffinity -8.5\n"
            "ATOM      1  C   LIG     1       0.000   0.000   0.000\n"
            "REMARK minimizedAffinity -7.2\n"
        )

        affinities = parse_smina_output(str(pdbqt))

        assert len(affinities) == 2
        assert affinities[0] == pytest.approx(-8.5)
        assert affinities[1] == pytest.approx(-7.2)


@pytest.mark.skipif(
    not AFFINITIES_CSV.exists(),
    reason="Standardized affinities CSV not available",
)
class TestGetTopHits:
    """Tests for get_top_hits()."""

    def test_returns_correct_number(self):
        df = get_top_hits(str(AFFINITIES_CSV), n=5)
        assert len(df) == 5

    def test_sorted_by_affinity(self):
        df = get_top_hits(str(AFFINITIES_CSV), n=10)
        affinities = df["affinity"].tolist()
        assert affinities == sorted(affinities)

    def test_standardized_top_five(self):
        df = get_top_hits(str(AFFINITIES_CSV), n=5)
        assert df["ligand_id"].tolist() == [
            "R11",
            "13U",
            "BAH",
            "12U",
            "607",
        ]


def test_standardized_parser_rejects_multiple_affinities(tmp_path):
    pdbqt = tmp_path / "multi_out.pdbqt"
    pdbqt.write_text(
        "REMARK minimizedAffinity -8.5\n"
        "ATOM      1  C   LIG     1       0.000   0.000   0.000\n"
        "REMARK minimizedAffinity -7.2\n"
    )
    assert parse_affinity(pdbqt) is None
