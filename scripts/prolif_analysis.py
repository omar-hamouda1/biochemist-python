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
ranking = ranking[
    ranking["status"].isin(["ok", "existing"])
].dropna(subset=["affinity"])
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

# ProLIF 2.2.x exposes the interaction metadata for a residue pair through
# Fingerprint.metadata(). Query every ligand/protein residue pair directly so
# valid interactions are not lost by the higher-level sparse IFP filtering.
fp = plf.Fingerprint()

print("\n2. Running ProLIF for each ligand...\n")

all_results = {}
ligand_residues = {}

for lid in LIGANDS:
    sdf_path = f"{FIXED_DIR}/{lid}_fixed.sd"

    if not os.path.exists(sdf_path):
        print(f"   SKIP {lid}: file not found")
        continue

    ligand_rdkit = Chem.MolFromMolFile(
        sdf_path,
        removeHs=False,
        sanitize=True,
    )

    if ligand_rdkit is None:
        print(f"   FAIL {lid}: cannot read SDF")
        continue

    ligand_mol = plf.Molecule.from_rdkit(ligand_rdkit)
    records = []
    residues = set()

    for ligand_resid, ligand_residue in ligand_mol.residues.items():
        for protein_resid, protein_residue in protein_mol.residues.items():
            metadata = fp.metadata(
                ligand_residue,
                protein_residue,
            )

            if not metadata:
                continue

            residues.add(protein_resid)

            for interaction_name, occurrences in metadata.items():
                # metadata values are tuples of occurrence dictionaries.
                for occurrence in occurrences:
                    records.append(
                        {
                            "ligand_residue": ligand_resid,
                            "protein_residue": protein_resid,
                            "interaction": interaction_name,
                            "metadata": occurrence,
                        }
                    )

    all_results[lid] = records
    ligand_residues[lid] = residues

    interaction_types = Counter(
        record["interaction"] for record in records
    )

    print(
        f"   {lid:5s}: {len(records)} interactions, "
        f"{len(residues)} residues"
    )

if not all_results:
    raise RuntimeError("No ProLIF results were generated")

print(f"\n{'=' * 70}")
print("  Detailed Interactions")
print(f"{'=' * 70}")

for lid in LIGANDS:
    if lid not in all_results:
        continue

    print(f"\n▶ {lid}:\n{'-' * 70}")

    records = all_results[lid]
    residues = ligand_residues[lid]
    interaction_types = Counter(
        record["interaction"] for record in records
    )

    print(f"   Protein Residues ({len(residues)}):")
    for residue in sorted(str(r) for r in residues):
        print(f"      - {residue}")

    print("   Interaction Types:")
    if interaction_types:
        for itype, count in interaction_types.most_common():
            print(f"      {itype:20s}: {count}")
    else:
        print("      None")

common = set.intersection(
    *(set(map(str, residues)) for residues in ligand_residues.values())
) if ligand_residues else set()

print(
    f"\n{'=' * 70}\n"
    "  Comparison: Common Residues Across Ligands\n"
    f"{'=' * 70}"
)

print(
    f"\n★ Common residues (found in ALL {len(ligand_residues)} ligands):"
)
for residue in sorted(common):
    print(f"   ✅ {residue}")

print("\n★ Residues per ligand:")
for lid in LIGANDS:
    residues = ligand_residues.get(lid, set())
    print(f"   {lid}: {sorted(map(str, residues))}")

summary_data = []

for lid in LIGANDS:
    records = all_results.get(lid, [])
    residues = ligand_residues.get(lid, set())
    interaction_types = Counter(
        record["interaction"] for record in records
    )

    summary_data.append(
        {
            "ligand_id": lid,
            "affinity": float(
                ranking.loc[
                    ranking["ligand_id"] == lid,
                    "affinity",
                ].iloc[0]
            ),
            "n_interactions": len(records),
            "n_residues": len(residues),
            "interaction_types": "; ".join(
                f"{name}:{count}"
                for name, count in interaction_types.most_common()
            ),
        }
    )

summary_df = (
    pd.DataFrame(summary_data)
    .sort_values("affinity")
    .reset_index(drop=True)
)

summary_df.to_csv(OUTPUT, index=False)

print(f"\n{summary_df.to_string(index=False)}")
print(f"\n   ✅ Saved: {OUTPUT}")
print(
    f"\n{'=' * 70}\n"
    "  🎉 Analysis Complete!\n"
    f"{'=' * 70}"
)
