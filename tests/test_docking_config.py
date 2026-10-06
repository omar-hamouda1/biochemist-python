"""Tests for the authoritative docking configuration."""

import os
from pathlib import Path

import pytest

from src.docking_config import DockingConfigError, load_docking_config
from scripts.run_standardized_docking import build_smina_command


CONFIG_PATH = Path("configs/docking_config.yml")


def test_loads_authoritative_config():
    config = load_docking_config(CONFIG_PATH)

    assert config.target_pdb == "2ZQ2"
    assert config.reference_ligand == "13U"
    assert config.center == (17.672, -8.256, 10.688)
    assert config.box_size == (25.0, 25.0, 25.0)
    assert config.exhaustiveness == 4
    assert config.num_modes == 1
    assert config.seed == 42
    assert config.timeout_seconds == 300
    assert config.expected_ligands == 111
    assert config.smina_executable == "smina"
    assert config.receptor == Path("docking/receptor/2zq2_receptor.pdbqt").resolve()


def test_default_config_is_independent_of_working_directory(tmp_path):
    original_cwd = Path.cwd()
    try:
        os.chdir(tmp_path)
        config = load_docking_config()
    finally:
        os.chdir(original_cwd)

    assert config.target_pdb == "2ZQ2"
    assert config.reference_ligand == "13U"
    assert config.receptor == Path("docking/receptor/2zq2_receptor.pdbqt").resolve()


def test_config_paths_are_project_root_relative(tmp_path):
    config_file = tmp_path / "configs" / "docking_config.yml"
    config_file.parent.mkdir()
    config_file.write_text(
        "\n".join(
            [
                "target_pdb: TEST",
                "reference_ligand: REF",
                "center_x: 1.0",
                "center_y: 2.0",
                "center_z: 3.0",
                "size_x: 10.0",
                "size_y: 11.0",
                "size_z: 12.0",
                "exhaustiveness: 2",
                "num_modes: 1",
                "seed: 7",
                "timeout_seconds: 60",
                "expected_ligands: 1",
                "receptor: docking/receptor.pdbqt",
                "ligand_dir: docking/ligands",
                "manifest: docking/manifest.csv",
                "output_dir: docking/results",
                "report: docking/report.csv",
            ]
        )
    )

    config = load_docking_config(config_file)
    assert config.receptor == tmp_path / "docking/receptor.pdbqt"
    assert config.output_dir == tmp_path / "docking/results"


def test_rejects_invalid_box_size(tmp_path):
    config_file = tmp_path / "configs" / "docking_config.yml"
    config_file.parent.mkdir()
    config_file.write_text(
        "\n".join(
            [
                "target_pdb: TEST",
                "reference_ligand: REF",
                "center_x: 1",
                "center_y: 2",
                "center_z: 3",
                "size_x: 0",
                "size_y: 10",
                "size_z: 10",
                "exhaustiveness: 2",
                "num_modes: 1",
                "seed: 7",
                "timeout_seconds: 60",
                "expected_ligands: 1",
                "receptor: receptor.pdbqt",
                "ligand_dir: ligands",
                "manifest: manifest.csv",
                "output_dir: results",
                "report: report.csv",
            ]
        )
    )

    with pytest.raises(DockingConfigError, match="size_x must be > 0"):
        load_docking_config(config_file)


def test_smina_command_uses_config_values():
    config = load_docking_config(CONFIG_PATH)
    ligand = config.ligand_dir / "R11.pdbqt"
    output = config.output_dir / "R11_out.pdbqt"

    command = build_smina_command(config, ligand, output)

    assert command[0] == "smina"
    assert "--center_x" in command
    assert command[command.index("--center_x") + 1] == "17.672"
    assert command[command.index("--exhaustiveness") + 1] == "4"
    assert command[command.index("--num_modes") + 1] == "1"
    assert command[command.index("--seed") + 1] == "42"


def test_rejects_non_integral_seed(tmp_path):
    config_file = tmp_path / "configs" / "docking_config.yml"
    config_file.parent.mkdir()
    config_file.write_text(
        "\n".join(
            [
                "target_pdb: TEST",
                "reference_ligand: REF",
                "center_x: 1",
                "center_y: 2",
                "center_z: 3",
                "size_x: 10",
                "size_y: 10",
                "size_z: 10",
                "exhaustiveness: 2",
                "num_modes: 1",
                "seed: 7.5",
                "timeout_seconds: 60",
                "expected_ligands: 1",
                "receptor: receptor.pdbqt",
                "ligand_dir: ligands",
                "manifest: manifest.csv",
                "output_dir: results",
                "report: report.csv",
            ]
        )
    )

    with pytest.raises(DockingConfigError, match="seed must be an integer"):
        load_docking_config(config_file)
