"""Prepare the validated ligand library for standardized docking."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import meeko
import rdkit
from rdkit import Chem

from src.docking_config import load_docking_config
from src.provenance import (
    build_artifact_metadata,
    executable_version,
    provenance_path,
    sha256_file,
    write_metadata,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
LIGANDS_DIR = PROJECT_ROOT / "ligands"
PREPARED_DIR = PROJECT_ROOT / "docking" / "ligands"
RESULTS_DIR = PROJECT_ROOT / "docking" / "results"
SCRIPT_PATH = Path(__file__).resolve()


def is_valid_pdbqt(path: Path) -> bool:
    """Return True only when a PDBQT exists and contains atoms."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open(encoding="utf-8", errors="replace") as handle:
        return any(line.startswith(("ATOM", "HETATM")) for line in handle)


def current_preparation_metadata(sdf_path: Path, method: str) -> dict:
    """Build provenance metadata for one prepared ligand."""
    config = load_docking_config()
    tool_versions = {
        "openbabel": executable_version("obabel"),
        "rdkit": rdkit.__version__,
        "meeko": meeko.__version__,
    }
    protocol = {
        "method": method,
        "hydrogens": True,
        "partial_charge": "gasteiger",
        "sdf_path": str(sdf_path.relative_to(PROJECT_ROOT)),
        "config_path": str(config.manifest.relative_to(PROJECT_ROOT)),
    }
    return build_artifact_metadata(
        source_path=sdf_path,
        producer_path=SCRIPT_PATH,
        tool_versions=tool_versions,
        protocol=protocol,
        artifact_kind="prepared_ligand_pdbqt",
    )


def valid_existing_preparation(sdf_path: Path, pdbqt_path: Path) -> bool:
    """Return True only when a prepared PDBQT has a matching provenance sidecar."""
    sidecar = provenance_path(pdbqt_path)
    if not is_valid_pdbqt(pdbqt_path) or not sidecar.exists():
        return False
    try:
        metadata = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False

    if metadata.get("source_sha256") != sha256_file(sdf_path):
        return False
    if metadata.get("producer_sha256") != sha256_file(SCRIPT_PATH):
        return False
    if metadata.get("output_sha256") != sha256_file(pdbqt_path):
        return False

    for method in ("openbabel", "meeko"):
        expected = current_preparation_metadata(sdf_path, method)
        if metadata.get("fingerprint") == expected.get("fingerprint"):
            return True
    return False


def write_preparation_sidecar(
    sdf_path: Path,
    pdbqt_path: Path,
    method: str,
) -> None:
    """Write provenance for a prepared PDBQT."""
    metadata = current_preparation_metadata(sdf_path, method)
    write_metadata(
        provenance_path(pdbqt_path),
        metadata,
        output_path=pdbqt_path,
    )


def run_openbabel(sdf_path: Path, output_path: Path):
    """Try SDF -> PDBQT using Open Babel."""
    if output_path.exists():
        output_path.unlink()

    command = [
        "obabel",
        "-isdf", str(sdf_path),
        "-opdbqt",
        "-O", str(output_path),
        "-h",
        "--partialcharge", "gasteiger",
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode == 0 and is_valid_pdbqt(output_path):
        return True, "Open Babel"

    message = result.stderr.strip() or result.stdout.strip()
    message = message.replace("\n", " ")[:200]
    return False, message or "Open Babel produced no valid PDBQT"


def run_meeko(sdf_path: Path, output_path: Path):
    """Try SDF -> PDBQT using Meeko."""
    if output_path.exists():
        output_path.unlink()

    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        if not supplier:
            return False, "RDKit could not read the SDF"

        mol = supplier[0]
        if mol is None:
            return False, "RDKit could not read the molecule"

        preparation = meeko.MoleculePreparation()
        setups = preparation.prepare(mol)
        if not setups:
            return False, "Meeko produced no setup"

        result = meeko.PDBQTWriterLegacy.write_string(setups[0])
        pdbqt = result[0] if isinstance(result, tuple) else result
        output_path.write_text(pdbqt, encoding="utf-8")

        if not is_valid_pdbqt(output_path):
            return False, "Meeko produced an invalid PDBQT"
        return True, "Meeko"
    except Exception as exc:
        return False, str(exc).replace("\n", " ")[:200]


def prepare_ligand(sdf_path: Path):
    """Prepare one ligand and retain it only with matching provenance."""
    output_path = PREPARED_DIR / f"{sdf_path.stem}.pdbqt"

    if valid_existing_preparation(sdf_path, output_path):
        return {
            "ligand_id": sdf_path.stem,
            "status": "already_valid",
            "method": "existing",
            "message": "Existing PDBQT kept with matching provenance",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    output_path.unlink(missing_ok=True)
    provenance_path(output_path).unlink(missing_ok=True)

    ok, message = run_openbabel(sdf_path, output_path)
    if ok:
        write_preparation_sidecar(sdf_path, output_path, "openbabel")
        return {
            "ligand_id": sdf_path.stem,
            "status": "prepared",
            "method": "openbabel",
            "message": "Open Babel succeeded",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    openbabel_error = message
    ok, message = run_meeko(sdf_path, output_path)
    if ok:
        write_preparation_sidecar(sdf_path, output_path, "meeko")
        return {
            "ligand_id": sdf_path.stem,
            "status": "prepared",
            "method": "meeko",
            "message": f"Open Babel failed: {openbabel_error}",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    return {
        "ligand_id": sdf_path.stem,
        "status": "failed",
        "method": "none",
        "message": f"Open Babel: {openbabel_error}; Meeko: {message}",
        "pdbqt_bytes": 0,
    }


def load_validated_sdf_files() -> list[Path]:
    """Return only SDF inputs listed in the authoritative validated manifest."""
    config = load_docking_config()
    with config.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    ligand_ids = [row["ligand_id"].strip() for row in rows if row.get("ligand_id")]
    if len(ligand_ids) != config.expected_ligands:
        raise RuntimeError(
            f"Expected {config.expected_ligands} validated ligands, found {len(ligand_ids)}"
        )
    if len(set(ligand_ids)) != len(ligand_ids):
        raise RuntimeError("Validated ligand manifest contains duplicate IDs")

    sdf_files = []
    for ligand_id in ligand_ids:
        sdf_path = LIGANDS_DIR / f"{ligand_id}.sdf"
        if not sdf_path.exists():
            raise FileNotFoundError(f"Validated ligand SDF not found: {sdf_path}")
        sdf_files.append(sdf_path)
    return sorted(sdf_files)


def main() -> None:
    """Prepare only the authoritative validated ligand set."""
    PREPARED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    sdf_files = load_validated_sdf_files()

    print("=" * 80)
    print("LIGAND PREPARATION")
    print("=" * 80)
    print(f"Total SDF ligands: {len(sdf_files)}")
    print()

    report = []
    for index, sdf_path in enumerate(sdf_files, start=1):
        result = prepare_ligand(sdf_path)
        report.append(result)
        print(
            f"[{index:3d}/{len(sdf_files)}] "
            f"{result['ligand_id']:>5} | "
            f"{result['status']:<13} | "
            f"{result['method']}"
        )

    report_path = RESULTS_DIR / "ligand_preparation.csv"
    with report_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand_id", "status", "method", "message", "pdbqt_bytes"],
        )
        writer.writeheader()
        writer.writerows(report)

    failed = sum(row["status"] == "failed" for row in report)
    print()
    print("=" * 80)
    print(f"Already valid : {sum(row['status'] == 'already_valid' for row in report)}")
    print(f"Prepared now  : {sum(row['status'] == 'prepared' for row in report)}")
    print(f"Failed        : {failed}")
    print(f"Report        : {report_path}")
    print("=" * 80)

    if failed:
        raise RuntimeError(f"Ligand preparation incomplete: {failed} failures")


if __name__ == "__main__":
    main()
