"""run_standardized_docking.py — Re-dock the validated ligand set under one fixed protocol.

This is the final homogeneous docking run used for scientific ranking.
The 111-ligand manifest comes from the recovered/validated dataset; the six
documented exceptions remain excluded from docking.
"""

from pathlib import Path
import csv
import subprocess


RECEPTOR = Path("docking/receptor/2zq2_receptor.pdbqt")
LIGAND_DIR = Path("docking/ligands")
RESULTS_DIR = Path("docking/results/standardized")
MANIFEST = Path("docking/results/affinities_recovered_111.csv")
REPORT = Path("docking/results/standardized_affinities.csv")

CENTER = (17.672, -8.256, 10.688)
BOX_SIZE = 25.0
EXHAUSTIVENESS = 4
NUM_MODES = 1
SEED = 42
TIMEOUT_SEC = 300


def parse_affinity(output_file: Path):
    """Read the best smina affinity from minimizedAffinity."""
    if not output_file.exists():
        return None

    with output_file.open() as handle:
        for line in handle:
            if line.startswith("REMARK minimizedAffinity"):
                try:
                    return float(line.split()[2])
                except (ValueError, IndexError):
                    return None
    return None


def load_ligands():
    """Load the 111-ligand validation manifest."""
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

    with MANIFEST.open(newline="") as handle:
        reader = csv.DictReader(handle)
        ligand_ids = [row["ligand_id"] for row in reader if row.get("ligand_id")]

    if len(ligand_ids) != 111:
        raise ValueError(
            f"Expected 111 validated ligands, found {len(ligand_ids)}"
        )

    if len(set(ligand_ids)) != len(ligand_ids):
        raise ValueError("Validated ligand manifest contains duplicate IDs")

    return ligand_ids


def run_one(ligand_id: str):
    """Run one standardized smina docking job."""
    ligand = LIGAND_DIR / f"{ligand_id}.pdbqt"
    output = RESULTS_DIR / f"{ligand_id}_out.pdbqt"

    if not ligand.exists():
        return {
            "ligand_id": ligand_id,
            "affinity": None,
            "status": "ligand_missing",
            "message": str(ligand),
        }

    if output.exists():
        output.unlink()

    command = [
        "smina",
        "--receptor",
        str(RECEPTOR),
        "--ligand",
        str(ligand),
        "--center_x",
        str(CENTER[0]),
        "--center_y",
        str(CENTER[1]),
        "--center_z",
        str(CENTER[2]),
        "--size_x",
        str(BOX_SIZE),
        "--size_y",
        str(BOX_SIZE),
        "--size_z",
        str(BOX_SIZE),
        "--exhaustiveness",
        str(EXHAUSTIVENESS),
        "--num_modes",
        str(NUM_MODES),
        "--seed",
        str(SEED),
        "--out",
        str(output),
    ]

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SEC,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "ligand_id": ligand_id,
            "affinity": None,
            "status": "timeout",
            "message": f">{TIMEOUT_SEC}s",
        }

    affinity = parse_affinity(output)

    if affinity is not None:
        return {
            "ligand_id": ligand_id,
            "affinity": affinity,
            "status": "ok",
            "message": "",
        }

    message = (result.stderr or result.stdout).replace("\n", " ").strip()
    return {
        "ligand_id": ligand_id,
        "affinity": None,
        "status": "failed",
        "message": message[:300],
    }


def main():
    if not RECEPTOR.exists():
        raise FileNotFoundError(f"Missing receptor: {RECEPTOR}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ligand_ids = load_ligands()

    print("=" * 80)
    print("STANDARDIZED DOCKING — 111 VALIDATED LIGANDS")
    print("=" * 80)
    print(f"Seed: {SEED}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"Num modes: {NUM_MODES}")
    print(f"Box: center={CENTER}, size={BOX_SIZE} Å")
    print(f"Ligands: {len(ligand_ids)}")
    print()

    results = []
    for index, ligand_id in enumerate(ligand_ids, start=1):
        row = run_one(ligand_id)
        results.append(row)
        print(
            f"[{index:3d}/{len(ligand_ids)}] "
            f"{ligand_id:>5} | {row['status']:<14} | "
            f"affinity={row['affinity']}"
        )

    with REPORT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ligand_id", "affinity", "status", "message"],
        )
        writer.writeheader()
        writer.writerows(results)

    successful = [r for r in results if r["status"] == "ok"]
    failed = [r for r in results if r["status"] != "ok"]

    print()
    print("=" * 80)
    print(f"Successful standardized dockings: {len(successful)}/{len(results)}")
    print(f"Failed/timeout: {len(failed)}")
    print(f"Report: {REPORT}")
    print(f"Pose directory: {RESULTS_DIR}")
    print("=" * 80)

    if failed:
        print("\nFailed ligands:")
        for row in failed:
            print(f"  {row['ligand_id']}: {row['status']} — {row['message']}")


if __name__ == "__main__":
    main()
