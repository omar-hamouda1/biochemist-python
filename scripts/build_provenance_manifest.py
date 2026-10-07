"""Build a tracked provenance manifest for the standardized docking run."""

from __future__ import annotations

import csv
import importlib.metadata
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from src.docking import parse_smina_output
from src.docking_config import load_docking_config
from src.provenance import (
    executable_version,
    load_metadata,
    provenance_path,
    sha256_file,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT_ROOT / "docking" / "results" / "standardized_provenance.json"
RUNNER = PROJECT_ROOT / "scripts" / "run_standardized_docking.py"
PARSER = PROJECT_ROOT / "src" / "docking.py"
PREPARER = PROJECT_ROOT / "scripts" / "prepare_ligands.py"
RECEPTOR_PREPARER = PROJECT_ROOT / "scripts" / "prepare_receptor.py"
PROTEIN_PREPARER = PROJECT_ROOT / "scripts" / "prepare_protein.py"


def package_version(name: str) -> str:
    """Return an installed package version without importing it."""
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def git_value(*args: str) -> str:
    """Return a Git value when the local checkout is a Git worktree."""
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"
    return result.stdout.strip()


def read_manifest_ids(path: Path) -> list[str]:
    """Read the validated ligand IDs."""
    with path.open(newline="", encoding="utf-8") as handle:
        return [row["ligand_id"] for row in csv.DictReader(handle)]


def main() -> None:
    """Validate current artifacts and write the tracked provenance manifest."""
    config = load_docking_config()

    if not config.receptor.exists():
        raise FileNotFoundError(f"Missing receptor: {config.receptor}")
    if not OUTPUT.parent.exists():
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    validated_ids = read_manifest_ids(config.manifest)
    report_rows = {
        row["ligand_id"]: row
        for row in csv.DictReader(config.report.open(newline="", encoding="utf-8"))
    }

    receptor_sidecar = provenance_path(config.receptor)
    receptor_metadata = load_metadata(receptor_sidecar)
    if receptor_metadata is None:
        raise RuntimeError(
            f"Missing or invalid receptor provenance: {receptor_sidecar}"
        )

    ligands = []
    docking_outputs = []
    errors = []

    for ligand_id in validated_ids:
        ligand = config.ligand_dir / f"{ligand_id}.pdbqt"
        ligand_sidecar = provenance_path(ligand)
        output = config.output_dir / f"{ligand_id}_out.pdbqt"
        output_sidecar = provenance_path(output)

        ligand_meta = load_metadata(ligand_sidecar)
        output_meta = load_metadata(output_sidecar)

        if ligand_meta is None:
            errors.append(f"{ligand_id}: missing ligand provenance")
        if output_meta is None:
            errors.append(f"{ligand_id}: missing docking provenance")
        if not ligand.exists():
            errors.append(f"{ligand_id}: missing prepared ligand")
        if not output.exists():
            errors.append(f"{ligand_id}: missing docking output")

        affinity_values = parse_smina_output(str(output))
        report = report_rows.get(ligand_id, {})
        if len(affinity_values) != 1:
            errors.append(
                f"{ligand_id}: expected exactly one affinity, found {len(affinity_values)}"
            )

        ligands.append(
            {
                "ligand_id": ligand_id,
                "source_sdf_sha256": sha256_file(
                    PROJECT_ROOT / "ligands" / f"{ligand_id}.sdf"
                ),
                "prepared_pdbqt_sha256": (
                    sha256_file(ligand) if ligand.exists() else None
                ),
                "preparation_fingerprint": (
                    ligand_meta.get("fingerprint") if ligand_meta else None
                ),
            }
        )
        docking_outputs.append(
            {
                "ligand_id": ligand_id,
                "affinity": (
                    float(report["affinity"])
                    if report.get("affinity") not in {None, ""}
                    else None
                ),
                "output_pdbqt_sha256": (
                    sha256_file(output) if output.exists() else None
                ),
                "docking_fingerprint": (
                    output_meta.get("fingerprint") if output_meta else None
                ),
                "recorded_output_sha256": (
                    output_meta.get("output_sha256") if output_meta else None
                ),
            }
        )

    if errors:
        raise RuntimeError(
            "Cannot build provenance manifest until all artifacts are valid:\n"
            + "\n".join(f"- {error}" for error in errors)
        )

    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit_at_generation": git_value("rev-parse", "HEAD"),
        "config": {
            "path": str(Path("configs/docking_config.yml")),
            "sha256": sha256_file(PROJECT_ROOT / "configs" / "docking_config.yml"),
        },
        "environment_snapshot": {
            "path": str(Path("docking/results/standardized_environment.txt")),
            "sha256": (
                sha256_file(PROJECT_ROOT / "docking" / "results" / "standardized_environment.txt")
                if (PROJECT_ROOT / "docking" / "results" / "standardized_environment.txt").exists()
                else None
            ),
        },
        "authoritative_files": {
            "validated_manifest": {
                "path": str(config.manifest.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(config.manifest),
            },
            "standardized_report": {
                "path": str(config.report.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(config.report),
            },
        },
        "receptor": {
            "path": str(config.receptor.relative_to(PROJECT_ROOT)),
            "sha256": sha256_file(config.receptor),
            "provenance": receptor_metadata,
            "preparation_script_sha256": sha256_file(RECEPTOR_PREPARER),
            "upstream_source_sha256": sha256_file(PROJECT_ROOT / "pdb" / "protein_h.pdb"),
            "upstream_preparation_script_sha256": sha256_file(PROTEIN_PREPARER),
            "upstream_notebook_sha256": sha256_file(
                PROJECT_ROOT / "notebooks" / "13_binding_site.ipynb"
            ),
        },
        "pipeline_code": {
            "runner_sha256": sha256_file(RUNNER),
            "parser_sha256": sha256_file(PARSER),
            "ligand_preparer_sha256": sha256_file(PREPARER),
            "receptor_preparer_sha256": sha256_file(RECEPTOR_PREPARER),
            "protein_preparer_sha256": sha256_file(PROTEIN_PREPARER),
        },
        "tools": {
            "smina": executable_version(config.smina_executable),
            "openbabel": executable_version("obabel"),
            "pdb2pqr": executable_version("pdb2pqr"),
            "python": sys.version.split()[0],
            "rdkit": package_version("rdkit"),
            "meeko": package_version("meeko"),
            "prolif": package_version("prolif"),
        },
        "validated_ligands": ligands,
        "docking_outputs": docking_outputs,
    }

    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote provenance manifest: {OUTPUT}")
    print(f"Validated ligands: {len(validated_ids)}")


if __name__ == "__main__":
    main()
