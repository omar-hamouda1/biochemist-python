from pathlib import Path
import csv
import subprocess

from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy


LIGANDS_DIR = Path("ligands")
PREPARED_DIR = Path("docking/ligands")
RESULTS_DIR = Path("docking/results")

PREPARED_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def is_valid_pdbqt(path: Path) -> bool:
    """Return True only when a PDBQT exists and contains atoms."""
    if not path.exists():
        return False

    if path.stat().st_size == 0:
        return False

    with path.open() as f:
        return any(
            line.startswith(("ATOM", "HETATM"))
            for line in f
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
    )

    if result.returncode == 0 and is_valid_pdbqt(output_path):
        return True, "Open Babel"

    message = result.stderr.strip() or result.stdout.strip()
    message = message.replace("\n", " ")[:200]

    if not message:
        message = "Open Babel produced no valid PDBQT"

    return False, message


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

        preparation = MoleculePreparation()
        setups = preparation.prepare(mol)

        if not setups:
            return False, "Meeko produced no setup"

        result = PDBQTWriterLegacy.write_string(setups[0])
        pdbqt = result[0] if isinstance(result, tuple) else result

        output_path.write_text(pdbqt)

        if not is_valid_pdbqt(output_path):
            return False, "Meeko produced an invalid PDBQT"

        return True, "Meeko"

    except Exception as exc:
        return False, str(exc).replace("\n", " ")[:200]


def prepare_ligand(sdf_path: Path):
    """Prepare one ligand, preserving already-valid PDBQT files."""
    ligand_id = sdf_path.stem
    output_path = PREPARED_DIR / f"{ligand_id}.pdbqt"

    if is_valid_pdbqt(output_path):
        return {
            "ligand_id": ligand_id,
            "status": "already_valid",
            "method": "existing",
            "message": "Existing PDBQT kept",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    if output_path.exists():
        output_path.unlink()

    ok, message = run_openbabel(sdf_path, output_path)

    if ok:
        return {
            "ligand_id": ligand_id,
            "status": "prepared",
            "method": "openbabel",
            "message": "Open Babel succeeded",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    openbabel_error = message

    ok, message = run_meeko(sdf_path, output_path)

    if ok:
        return {
            "ligand_id": ligand_id,
            "status": "prepared",
            "method": "meeko",
            "message": f"Open Babel failed: {openbabel_error}",
            "pdbqt_bytes": output_path.stat().st_size,
        }

    return {
        "ligand_id": ligand_id,
        "status": "failed",
        "method": "none",
        "message": (
            f"Open Babel: {openbabel_error}; "
            f"Meeko: {message}"
        ),
        "pdbqt_bytes": 0,
    }


def main():
    sdf_files = sorted(LIGANDS_DIR.glob("*.sdf"))

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

    with report_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "ligand_id",
                "status",
                "method",
                "message",
                "pdbqt_bytes",
            ],
        )
        writer.writeheader()
        writer.writerows(report)

    prepared = sum(r["status"] == "prepared" for r in report)
    existing = sum(r["status"] == "already_valid" for r in report)
    failed = sum(r["status"] == "failed" for r in report)

    print()
    print("=" * 80)
    print(f"Already valid : {existing}")
    print(f"Prepared now  : {prepared}")
    print(f"Failed        : {failed}")
    print(f"Report        : {report_path}")
    print("=" * 80)


if __name__ == "__main__":
    main()
