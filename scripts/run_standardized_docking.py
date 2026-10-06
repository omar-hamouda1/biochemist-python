"""Resumable standardized docking for the validated ligand set."""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

from src.docking import parse_smina_output
from src.docking_config import DockingConfig, load_docking_config


def build_smina_command(
    config: DockingConfig,
    ligand: Path,
    output: Path,
) -> list[str]:
    """Build the Smina command directly from the validated config."""

    return [
        config.smina_executable,
        "--receptor", str(config.receptor),
        "--ligand", str(ligand),
        "--center_x", str(config.center_x),
        "--center_y", str(config.center_y),
        "--center_z", str(config.center_z),
        "--size_x", str(config.size_x),
        "--size_y", str(config.size_y),
        "--size_z", str(config.size_z),
        "--exhaustiveness", str(config.exhaustiveness),
        "--num_modes", str(config.num_modes),
        "--seed", str(config.seed),
        "--out", str(output),
    ]


def parse_affinity(path: Path) -> float | None:
    """Return the first affinity recorded in a Smina output file."""

    affinities = parse_smina_output(str(path))
    return affinities[0] if affinities else None


def load_ligands(config: DockingConfig) -> list[str]:
    """Load and validate ligand IDs from the authoritative manifest."""

    with config.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    ligand_ids = [row["ligand_id"] for row in rows if row.get("ligand_id")]
    if len(ligand_ids) != config.expected_ligands:
        raise ValueError(
            f"Expected {config.expected_ligands} validated ligands, "
            f"found {len(ligand_ids)}"
        )
    if len(set(ligand_ids)) != config.expected_ligands:
        raise ValueError("Validated ligand manifest contains duplicate IDs")
    return ligand_ids


def run_one(config: DockingConfig, ligand_id: str):
    ligand = config.ligand_dir / f"{ligand_id}.pdbqt"
    output = config.output_dir / f"{ligand_id}_out.pdbqt"

    if not ligand.exists():
        return ligand_id, None, "ligand_missing", str(ligand)

    existing = parse_affinity(output)
    if existing is not None:
        return ligand_id, existing, "existing", ""

    if output.exists():
        output.unlink()

    command = build_smina_command(config, ligand, output)

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=config.timeout_seconds,
            check=False,
        )
    except FileNotFoundError:
        return ligand_id, None, "smina_missing", config.smina_executable
    except subprocess.TimeoutExpired:
        return ligand_id, None, "timeout", f">{config.timeout_seconds}s"

    affinity = parse_affinity(output)
    if affinity is not None:
        return ligand_id, affinity, "ok", ""

    message = (result.stderr or result.stdout).replace("\n", " ").strip()
    return ligand_id, None, "failed", message[:300]


def main() -> None:
    config = load_docking_config()

    if not config.receptor.exists():
        raise FileNotFoundError(f"Missing receptor: {config.receptor}")
    if shutil.which(config.smina_executable) is None:
        raise FileNotFoundError(
            f"Smina executable not found on PATH: {config.smina_executable}"
        )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    ligand_ids = load_ligands(config)
    rows = []

    print("=" * 80)
    print(f"STANDARDIZED DOCKING — RESUMABLE {config.expected_ligands}-LIGAND RUN")
    print("=" * 80)

    for index, ligand_id in enumerate(ligand_ids, start=1):
        lid, affinity, status, message = run_one(config, ligand_id)
        rows.append(
            {
                "ligand_id": lid,
                "affinity": affinity,
                "status": status,
                "message": message,
            }
        )
        print(
            f"[{index:3d}/{len(ligand_ids)}] {lid:>5} | "
            f"{status:<14} | affinity={affinity}",
            flush=True,
        )

    with config.report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand_id", "affinity", "status", "message"],
        )
        writer.writeheader()
        writer.writerows(rows)

    completed = [r for r in rows if r["status"] in {"ok", "existing"}]
    newly = [r for r in rows if r["status"] == "ok"]
    reused = [r for r in rows if r["status"] == "existing"]
    failed = [r for r in rows if r["status"] not in {"ok", "existing"}]

    print("=" * 80)
    print(f"Completed standardized dockings: {len(completed)}/{len(rows)}")
    print(f"Newly docked this run: {len(newly)}")
    print(f"Reused existing valid outputs: {len(reused)}")
    print(f"Failed/timeout: {len(failed)}")
    print(f"Report: {config.report}")
    print(f"Pose directory: {config.output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
