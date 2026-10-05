import os
from collections import Counter

import MDAnalysis as mda
import pandas as pd
import prolif as plf
from rdkit import Chem

DOCKING_RESULTS = "docking/results/standardized_affinities.csv"
PROTEIN = "pdb/protein_h.pdb"
FIXED_DIR = "docking/results/fixed"
OUTPUT = "docking/results/prolif_summary.csv"

ranking = pd.read_csv(DOCKING_RESULTS)
ranking["affinity"] = pd.to_numeric(ranking["affinity"], errors="coerce")
ranking = ranking[ranking["status"].isin(["ok", "existing"])].dropna(subset=["affinity"])
ranking = ranking.sort_values("affinity").reset_index(drop=True)
LIGANDS = ranking.head(5)["ligand_id"].tolist()

print("=" * 70)
print("  ProLIF Analysis: Standardized Top 5 Docked Ligands")
print("=" * 70)
print(f"\nTop 5: {LIGANDS}")

print("\n1. Loading protein...")
u = mda.Universe(PROTEIN)
protein_atoms = u.select_atoms("protein")
protein_mol = plf.Molecule.from_mda(protein_atoms)
print(f"   Protein: {len(protein_atoms)} atoms")

print("\n2. Running ProLIF for each ligand...\n")
all_results = {}
for lid in LIGANDS:
    sdf_path = f"{FIXED_DIR}/{lid}_fixed.sd"
    if not os.path.exists(sdf_path):
        print(f"   SKIP {lid}: file not found")
        continue
    ligand_rdkit = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if ligand_rdkit is None:
        print(f"   FAIL {lid}: cannot read SDF")
        continue
    ligand_mol = plf.Molecule.from_rdkit(ligand_rdkit)
    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand_mol], protein_mol)
    df = fp.to_dataframe()
    all_results[lid] = df
    print(f"   {lid:5s}: {df.shape[1]} interactions")

if not all_results:
    raise RuntimeError("No ProLIF results were generated")

print(f"\n{'=' * 70}")
print("  Detailed Interactions")
print(f"{'=' * 70}")

ligand_residues = {}
for lid, df in all_results.items():
    print(f"\n▶ {lid}:\n{'-' * 70}")
    residues = set()
    interaction_types = Counter()
    for col in df.columns:
        if len(col) >= 3:
            protein_res, itype = col[1], col[2]
        elif len(col) == 2:
            protein_res, itype = col[1], "Interaction"
        else:
            continue
        if df.iloc[0][col]:
            residues.add(protein_res)
            interaction_types[itype] += 1
    ligand_residues[lid] = residues
    print(f"   Protein Residues ({len(residues)}):")
    for r in sorted(residues):
        print(f"      - {r}")
    print("   Interaction Types:")
    for itype, count in interaction_types.most_common():
        print(f"      {itype:20s}: {count}")

common = set.intersection(*ligand_residues.values()) if ligand_residues else set()
print(f"\n{'=' * 70}\n  Comparison: Common Residues Across Ligands\n{'=' * 70}")
print(f"\n★ Common residues (found in ALL {len(ligand_residues)} ligands):")
for r in sorted(common):
    print(f"   ✅ {r}")
print("\n★ Residues per ligand:")
for lid, residues in ligand_residues.items():
    print(f"   {lid}: {sorted(residues)}")

summary_data = []
for lid, df in all_results.items():
    summary_data.append({
        "ligand_id": lid,
        "affinity": float(ranking.loc[ranking["ligand_id"] == lid, "affinity"].iloc[0]),
        "n_interactions": df.shape[1],
        "n_residues": len(ligand_residues[lid]),
    })
summary_df = pd.DataFrame(summary_data).sort_values("affinity")
summary_df.to_csv(OUTPUT, index=False)
print(f"\n{summary_df.to_string(index=False)}")
print(f"\n   ✅ Saved: {OUTPUT}")
print(f"\n{'=' * 70}\n  🎉 Analysis Complete!\n{'=' * 70}")
