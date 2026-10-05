"""Resumable standardized docking for the validated 111-ligand set."""

from pathlib import Path
import csv
import subprocess

RECEPTOR = Path("docking/receptor/2zq2_receptor.pdbqt")
LIGAND_DIR = Path("docking/ligands")
RESULTS_DIR = Path("docking/results/standardized")
MANIFEST = Path("docking/results/validated_ligands_111.csv")
REPORT = Path("docking/results/standardized_affinities.csv")

CENTER = (17.672, -8.256, 10.688)
BOX_SIZE = 25.0
EXHAUSTIVENESS = 4
NUM_MODES = 1
SEED = 42
TIMEOUT_SEC = 300


def parse_affinity(path: Path):
    if not path.exists():
        return None
    with path.open() as handle:
        for line in handle:
            if line.startswith("REMARK minimizedAffinity"):
                try:
                    return float(line.split()[2])
                except (ValueError, IndexError):
                    return None
    return None


def load_ligands():
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ligand_ids = [row["ligand_id"] for row in rows if row.get("ligand_id")]
    if len(ligand_ids) != 111:
        raise ValueError(f"Expected 111 validated ligands, found {len(ligand_ids)}")
    if len(set(ligand_ids)) != 111:
        raise ValueError("Validated ligand manifest contains duplicate IDs")
    return ligand_ids


def run_one(ligand_id: str):
    ligand = LIGAND_DIR / f"{ligand_id}.pdbqt"
    output = RESULTS_DIR / f"{ligand_id}_out.pdbqt"

    if not ligand.exists():
        return ligand_id, None, "ligand_missing", str(ligand)

    existing = parse_affinity(output)
    if existing is not None:
        return ligand_id, existing, "existing", ""

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

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_SEC,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return ligand_id, None, "timeout", f">{TIMEOUT_SEC}s"

    affinity = parse_affinity(output)
    if affinity is not None:
        return ligand_id, affinity, "ok", ""

    message = (result.stderr or result.stdout).replace("\n", " ").strip()
    return ligand_id, None, "failed", message[:300]


def main():
    if not RECEPTOR.exists():
        raise FileNotFoundError(f"Missing receptor: {RECEPTOR}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ligand_ids = load_ligands()
    rows = []

    print("=" * 80)
    print("STANDARDIZED DOCKING — RESUMABLE 111-LIGAND RUN")
    print("=" * 80)

    for index, ligand_id in enumerate(ligand_ids, start=1):
        lid, affinity, status, message = run_one(ligand_id)
        rows.append({
            "ligand_id": lid,
            "affinity": affinity,
            "status": status,
            "message": message,
        })
        print(
            f"[{index:3d}/{len(ligand_ids)}] {lid:>5} | "
            f"{status:<14} | affinity={affinity}",
            flush=True,
        )

    with REPORT.open("w", newline="") as handle:
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
    print(f"Report: {REPORT}")
    print(f"Pose directory: {RESULTS_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    main()