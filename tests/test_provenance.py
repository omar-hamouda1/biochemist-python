"""Tests for docking provenance metadata."""

import json
from pathlib import Path

from src.provenance import (
    build_fingerprint,
    build_metadata,
    load_metadata,
    provenance_path,
    sha256_file,
    validate_existing_output,
    write_metadata,
)


def _metadata(tmp_path):
    config = tmp_path / "config.yml"
    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    runner = tmp_path / "runner.py"
    parser = tmp_path / "parser.py"
    for path, content in (
        (config, "seed: 42\n"),
        (receptor, "ATOM receptor\n"),
        (ligand, "ATOM ligand\n"),
        (runner, "print('runner')\n"),
        (parser, "print('parser')\n"),
    ):
        path.write_text(content)

    command = ["smina", "--seed", "42"]
    metadata = build_metadata(
        config_path=config,
        receptor_path=receptor,
        ligand_path=ligand,
        runner_path=runner,
        parser_path=parser,
        smina_version="smina 2020.12.10",
        command=command,
    )
    return metadata


def test_fingerprint_is_stable(tmp_path):
    metadata = _metadata(tmp_path)
    assert metadata["fingerprint"] == build_fingerprint(metadata)


def test_sidecar_round_trip_and_output_hash(tmp_path):
    output = tmp_path / "pose.pdbqt"
    output.write_text("ATOM pose\n")
    sidecar = provenance_path(output)
    metadata = _metadata(tmp_path)

    write_metadata(
        sidecar,
        metadata,
        output_sha256=sha256_file(output),
        affinity=-8.5,
    )

    loaded = load_metadata(sidecar)
    assert loaded is not None
    assert loaded["output_sha256"] == sha256_file(output)
    assert loaded["affinity"] == -8.5
    assert validate_existing_output(
        output_path=output,
        metadata_path=sidecar,
        expected_metadata=metadata,
        expected_affinity=-8.5,
    )


def test_existing_output_is_rejected_when_input_fingerprint_changes(tmp_path):
    output = tmp_path / "pose.pdbqt"
    output.write_text("ATOM pose\n")
    sidecar = provenance_path(output)
    metadata = _metadata(tmp_path)
    write_metadata(
        sidecar,
        metadata,
        output_sha256=sha256_file(output),
        affinity=-8.5,
    )

    changed_config = tmp_path / "config.yml"
    changed_config.write_text("seed: 43\n")
    changed = dict(metadata)
    changed["config_sha256"] = sha256_file(changed_config)
    changed["fingerprint"] = build_fingerprint(changed)

    assert not validate_existing_output(
        output_path=output,
        metadata_path=sidecar,
        expected_metadata=changed,
        expected_affinity=-8.5,
    )


def test_existing_output_is_rejected_when_pose_changes(tmp_path):
    output = tmp_path / "pose.pdbqt"
    output.write_text("ATOM pose\n")
    sidecar = provenance_path(output)
    metadata = _metadata(tmp_path)
    write_metadata(
        sidecar,
        metadata,
        output_sha256=sha256_file(output),
        affinity=-8.5,
    )

    output.write_text("ATOM changed\n")

    assert not validate_existing_output(
        output_path=output,
        metadata_path=sidecar,
        expected_metadata=metadata,
        expected_affinity=-8.5,
    )


def test_fingerprint_ignores_machine_specific_command_paths(tmp_path):
    metadata = _metadata(tmp_path)
    alternate = dict(metadata)
    alternate["protocol"] = dict(metadata["protocol"])
    alternate["protocol"]["command"] = [
        "/different/conda/env/bin/smina",
        "--receptor",
        "/other/worktree/receptor.pdbqt",
    ]
    alternate["fingerprint"] = build_fingerprint(alternate)

    assert alternate["fingerprint"] == metadata["fingerprint"]
