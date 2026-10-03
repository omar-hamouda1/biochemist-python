import os
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
from collections import Counter
import pandas as pd

print("=" * 70)
print("  ProLIF Analysis: Top 4 Docked Ligands")
print("=" * 70)

PROTEIN = "pdb/protein_h.pdb"
FIXED_DIR = "docking/results/fixed"
LIGANDS = ["13U", "BAH", "607", "12U"]

# ─── Load protein ───
print("\n1. Loading protein...")
u = mda.Universe(PROTEIN)
protein_atoms = u.select_atoms("protein")
protein_mol = plf.Molecule.from_mda(protein_atoms)
print(f"   Protein: {len(protein_atoms)} atoms")

# ─── Run ProLIF for each ligand ───
print("\n2. Running ProLIF for each ligand...\n")

all_results = {}

for lid in LIGANDS:
    sdf_path = f"{FIXED_DIR}/{lid}_fixed.sdf"
    
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

# ─── Detailed Analysis (FIXED) ───
print(f"\n{'='*70}")
print(f"  Detailed Interactions")
print(f"{'='*70}")

for lid, df in all_results.items():
    print(f"\n▶ {lid}:")
    print("-" * 70)
    
    # Get protein residues (col[1] is the protein residue)
    residues = set()
    interaction_types = Counter()
    
    for col in df.columns:
        # ProLIF 2.x format: (ligand_res, protein_res, interaction)
        if len(col) >= 3:
            protein_res = col[1]
            itype = col[2]
        elif len(col) == 2:
            protein_res = col[1]
            itype = "Interaction"
        else:
            continue
        
        if df.iloc[0][col]:
            residues.add(protein_res)
            interaction_types[itype] += 1
    
    print(f"   Protein Residues ({len(residues)}):")
    for r in sorted(residues):
        print(f"      - {r}")
    
    print(f"   Interaction Types:")
    for itype, count in interaction_types.most_common():
        print(f"      {itype:20s}: {count}")

# ─── Compare all ligands ───
print(f"\n{'='*70}")
print(f"  Comparison: Common Residues Across Ligands")
print(f"{'='*70}")

ligand_residues = {}
for lid, df in all_results.items():
    residues = set()
    for col in df.columns:
        if len(col) >= 2:
            protein_res = col[1]
            if df.iloc[0][col]:
                residues.add(protein_res)
    ligand_residues[lid] = residues

common = set.intersection(*ligand_residues.values()) if ligand_residues else set()

print(f"\n★ Common residues (found in ALL {len(ligand_residues)} ligands):")
for r in sorted(common):
    print(f"   ✅ {r}")

print(f"\n★ Residues per ligand:")
for lid, residues in ligand_residues.items():
    print(f"   {lid}: {sorted(residues)}")

# ─── Summary Table ───
print(f"\n{'='*70}")
print(f"  Summary")
print(f"{'='*70}")

summary_data = []
for lid, df in all_results.items():
    residues = ligand_residues[lid]
    
    aff_df = pd.read_csv("docking/results/affinities_full.csv")
    aff = aff_df[aff_df["ligand_id"] == lid]["affinity"].values
    affinity = aff[0] if len(aff) > 0 else None
    
    summary_data.append({
        "ligand_id": lid,
        "affinity": affinity,
        "n_interactions": df.shape[1],
        "n_residues": len(residues),
    })

summary_df = pd.DataFrame(summary_data)
summary_df = summary_df.sort_values("affinity")
summary_path = "docking/results/prolif_summary.csv"
summary_df.to_csv(summary_path, index=False)

print(f"\n{summary_df.to_string(index=False)}")
print(f"\n   ✅ Saved: {summary_path}")

print(f"\n{'='*70}")
print(f"  🎉 Analysis Complete!")
print(f"{'='*70}")