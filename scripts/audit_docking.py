from pathlib import Path


LIGANDS_DIR = Path("ligands")
PREPARED_DIR = Path("docking/ligands")
RESULTS_DIR = Path("docking/results")


sdf_files = sorted(LIGANDS_DIR.glob("*.sdf"))

print(f"Total ligands found: {len(sdf_files)}")


preparation_failed = []
docking_failed = []


for sdf in sdf_files:
    ligand_id = sdf.stem

    pdbqt = PREPARED_DIR / f"{ligand_id}.pdbqt"
    output = RESULTS_DIR / f"{ligand_id}_out.pdbqt"

    preparation_ok = pdbqt.exists() and pdbqt.stat().st_size > 0
    docking_ok = output.exists() and output.stat().st_size > 0

    if not preparation_ok:
        preparation_failed.append(ligand_id)

    if not docking_ok:
        docking_failed.append(ligand_id)


print(f"Preparation failed: {len(preparation_failed)}")
print(f"Docking output failed: {len(docking_failed)}")

print("\nPreparation failures:")
for ligand_id in preparation_failed:
    print(f"  {ligand_id}")

print("\nDocking output failures:")
for ligand_id in docking_failed:
    print(f"  {ligand_id}")


import csv


affinities_file = RESULTS_DIR / "affinities_full.csv"

result_ids = set()

if affinities_file.exists():
    with affinities_file.open(newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            ligand_id = row.get("ligand_id")

            if ligand_id:
                result_ids.add(ligand_id)


sdf_ids = {sdf.stem for sdf in sdf_files}

missing_from_results = sorted(sdf_ids - result_ids)


print("\nResults CSV audit:")
print(f"  Ligands in SDF: {len(sdf_ids)}")
print(f"  Ligands in affinities_full.csv: {len(result_ids)}")
print(f"  Missing from results CSV: {len(missing_from_results)}")

for ligand_id in missing_from_results:
    print(f"  {ligand_id}")
