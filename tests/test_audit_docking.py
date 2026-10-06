"""Tests for the authoritative standardized docking audit."""

import csv

from src.docking_config import load_docking_config
from scripts.audit_docking import audit_standardized_docking


def _write_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _make_project(tmp_path, *, report_rows=None, exception_rows=None):
    (tmp_path / "configs").mkdir()
    (tmp_path / "ligands").mkdir()
    (tmp_path / "docking" / "results").mkdir(parents=True)
    (tmp_path / "docking" / "ligands").mkdir(parents=True)
    (tmp_path / "docking" / "results" / "standardized").mkdir(parents=True)

    (tmp_path / "configs" / "docking_config.yml").write_text(
        "\n".join(
            [
                "target_pdb: 2ZQ2",
                "reference_ligand: 13U",
                "center_x: 1",
                "center_y: 2",
                "center_z: 3",
                "size_x: 10",
                "size_y: 10",
                "size_z: 10",
                "smina_executable: smina",
                "exhaustiveness: 4",
                "num_modes: 1",
                "seed: 42",
                "timeout_seconds: 300",
                "expected_ligands: 2",
                "receptor: docking/receptor.pdbqt",
                "ligand_dir: docking/ligands",
                "manifest: docking/results/validated.csv",
                "output_dir: docking/results/standardized",
                "report: docking/results/standardized_affinities.csv",
            ]
        )
    )

    for ligand_id in ("AAA", "BBB", "CCC"):
        (tmp_path / "ligands" / f"{ligand_id}.sdf").write_text("placeholder\n")

    _write_csv(
        tmp_path / "docking" / "results" / "validated.csv",
        ["ligand_id"],
        [{"ligand_id": "AAA"}, {"ligand_id": "BBB"}],
    )
    _write_csv(
        tmp_path / "docking" / "results" / "docking_exceptions.csv",
        ["ligand_id", "stage", "status", "reason"],
        exception_rows
        or [
            {
                "ligand_id": "CCC",
                "stage": "preparation",
                "status": "unresolved",
                "reason": "test",
            }
        ],
    )
    _write_csv(
        tmp_path / "docking" / "results" / "standardized_affinities.csv",
        ["ligand_id", "affinity", "status", "message"],
        report_rows
        or [
            {"ligand_id": "AAA", "affinity": "-8.1", "status": "ok", "message": ""},
            {"ligand_id": "BBB", "affinity": "-7.9", "status": "ok", "message": ""},
        ],
    )

    return load_docking_config(tmp_path / "configs" / "docking_config.yml")


def test_authoritative_audit_passes(tmp_path):
    config = _make_project(tmp_path)
    audit = audit_standardized_docking(config)

    assert audit.ok
    assert len(audit.candidate_ids) == 3
    assert len(audit.validated_ids) == 2
    assert audit.exception_ids == ("CCC",)
    assert audit.report_ids == ("AAA", "BBB")


def test_audit_rejects_report_missing_validated_ligand(tmp_path):
    config = _make_project(
        tmp_path,
        report_rows=[
            {"ligand_id": "AAA", "affinity": "-8.1", "status": "ok", "message": ""},
        ],
    )
    audit = audit_standardized_docking(config)

    assert not audit.ok
    assert any("missing from standardized report" in error for error in audit.errors)


def test_audit_rejects_unknown_report_ligand(tmp_path):
    config = _make_project(
        tmp_path,
        report_rows=[
            {"ligand_id": "AAA", "affinity": "-8.1", "status": "ok", "message": ""},
            {"ligand_id": "ZZZ", "affinity": "-7.9", "status": "ok", "message": ""},
        ],
    )
    audit = audit_standardized_docking(config)

    assert not audit.ok
    assert any("Unknown ligands in standardized report" in error for error in audit.errors)


def test_audit_rejects_overlap_between_manifest_and_exceptions(tmp_path):
    config = _make_project(
        tmp_path,
        exception_rows=[
            {
                "ligand_id": "AAA",
                "stage": "preparation",
                "status": "unresolved",
                "reason": "test",
            }
        ],
    )
    audit = audit_standardized_docking(config)

    assert not audit.ok
    assert any("both validated manifest and exceptions" in error for error in audit.errors)


def test_audit_rejects_non_complete_report_status(tmp_path):
    config = _make_project(
        tmp_path,
        report_rows=[
            {"ligand_id": "AAA", "affinity": "-8.1", "status": "failed", "message": "x"},
            {"ligand_id": "BBB", "affinity": "-7.9", "status": "ok", "message": ""},
        ],
    )
    audit = audit_standardized_docking(config)

    assert not audit.ok
    assert any("non-complete status" in error for error in audit.errors)


def test_audit_rejects_non_numeric_affinity(tmp_path):
    config = _make_project(
        tmp_path,
        report_rows=[
            {"ligand_id": "AAA", "affinity": "bad", "status": "ok", "message": ""},
            {"ligand_id": "BBB", "affinity": "-7.9", "status": "ok", "message": ""},
        ],
    )
    audit = audit_standardized_docking(config)

    assert not audit.ok
    assert any("Non-numeric affinity" in error for error in audit.errors)
