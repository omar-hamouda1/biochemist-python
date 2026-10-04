from pathlib import Path
import csv
import subprocess


RECEPTOR = Path("docking/receptor/2zq2_receptor.pdbqt")
RESULTS_DIR = Path("docking/results")

CENTER = (17.672, -8.256, 10.688)
BOX_SIZE = 25.0

EXHAUSTIVENESS = 4
NUM_MODES = 1
SEED = 42

LIGANDS = [
    "132", "334", "678", "762", "847",
    "907", "972", "991", "BR6", "BRV", "ZEN",
]


def parse_affinity(output_file: Path):
    """Read minimizedAffinity from a smina PDBQT output file."""
    if not output_file.exists():
        return None

    with output_file.open() as f:
        for line in f:
            if line.startswith("REMARK minimizedAffinity"):
                try:
                    return float(line.split()[2])
                except (ValueError, IndexError):
                    return None

    return None


def run_docking(ligand_id: str):
    """Run one reproducible recovery docking."""
    ligand = Path(f"docking/ligands/{ligand_id}.pdbqt")
    output = RESULTS_DIR / f"{ligand_id}_out.pdbqt"

    if output.exists():
        output.unlink()

    command = [
        "smina",
        "--receptor", str(RECEPTOR),
        "--ligand", str(ligand),
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(BOX_SIZE),
        "--size_y", str(BOX_SIZE),
        "--size_z", str(BOX_SIZE),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--seed", str(SEED),
        "--out", str(output),
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=300,
    )

    affinity = parse_affinity(output)

    if affinity is not None:
        return {
            "ligand_id": ligand_id,
            "affinity": affinity,
            "status": "ok",
        }

    message = (result.stderr or result.stdout).strip()
    message = message.replace("\n", " ")[:300]

    return {
        "ligand_id": ligand_id,
        "affinity": None,
        "status": "failed",
        "message": message,
    }


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("RECOVERY DOCKING")
    print("=" * 80)
    print(f"Seed: {SEED}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print(f"Ligands: {len(LIGANDS)}")
    print()

    results = []

    for i, ligand_id in enumerate(LIGANDS, start=1):
        result = run_docking(ligand_id)
        results.append(result)

        print(
            f"[{i:2d}/{len(LIGANDS)}] "
            f"{ligand_id:>5} | "
            f"{result['status']} | "
            f"affinity={result['affinity']}"
        )

    report = RESULTS_DIR / "recovered_affinities.csv"

    with report.open("w", newline="") as f:
        fieldnames = ["ligand_id", "affinity", "status"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    successful = sum(r["status"] == "ok" for r in results)

    print()
    print("=" * 80)
    print(f"Successful recoveries: {successful}/{len(results)}")
    print(f"Report: {report}")
    print("=" * 80)


if __name__ == "__main__":
    main()
