"""Run ProLIF analysis for the current standardized docking Top 5."""

from collections import Counter
from pathlib import Path

import MDAnalysis as mda
import pandas as pd
import prolif as plf
from rdkit import Chem

from src.docking_config import load_docking_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TOP_N = 5


def load_top_hits(report_path: Path):
    ranking = pd.read_csv(report_path)
    required = {"ligand_id", "affinity", "status"}
    missing = required - set(ranking.columns)
    if missing:
        raise RuntimeError(
            f"Standardized report missing columns: {sorted(missing)}"
        )

    ranking["affinity"] = pd.to_numeric(
        ranking["affinity"], errors="coerce"
    )
    ranking = ranking[
        ranking["status"].isin(["ok", "existing"])
    ].dropna(subset=["affinity"])
    ranking = ranking.sort_values("affinity").reset_index(drop=True)

    top = ranking.head(TOP_N)
    if len(top) != TOP_N or top["ligand_id"].duplicated().any():
        raise RuntimeError("Current standardized Top 5 is incomplete or duplicated")

    return top


def main():
    config = load_docking_config()
    ranking = load_top_hits(config.report)
    ligands = ranking["ligand_id"].astype(str).tolist()

    protein_path = PROJECT_ROOT / "pdb" / "protein_h.pdb"
    fixed_dir = PROJECT_ROOT / "docking" / "results" / "fixed"
    output = PROJECT_ROOT / "docking" / "results" / "prolif_summary.csv"

    if not protein_path.exists():
        raise FileNotFoundError(f"Missing prepared protein: {protein_path}")

    missing_fixed = [
        lid for lid in ligands
        if not (fixed_dir / f"{lid}_fixed.sd").exists()
    ]
    if missing_fixed:
        raise RuntimeError(
            f"Missing ProLIF-ready Top 5 poses: {missing_fixed}. "
            "Run scripts/prepare_top5.py first."
        )

    print("=" * 70)
    print("  ProLIF Analysis: Standardized Top 5")
    print("=" * 70)
    print(f"\nTop 5: {ligands}")

    u = mda.Universe(str(protein_path))
    protein_atoms = u.select_atoms("protein")
    protein_mol = plf.Molecule.from_mda(protein_atoms)

    fp = plf.Fingerprint()
    protein_ids = [str(resid) for resid in protein_mol.residues.keys()]

    all_results = {}
    ligand_residues = {}

    for lid in ligands:
        sdf_path = fixed_dir / f"{lid}_fixed.sd"
        ligand_rdkit = Chem.MolFromMolFile(
            str(sdf_path), removeHs=False, sanitize=True
        )

        if ligand_rdkit is None:
            raise RuntimeError(f"RDKit could not parse {sdf_path}")

        ligand_mol = plf.Molecule.from_rdkit(ligand_rdkit)

        if ligand_mol.n_residues != 1:
            raise RuntimeError(
                f"{lid} contains {ligand_mol.n_residues} ligand residues; "
                "expected exactly one."
            )

        ligand_residue = ligand_mol[0]
        records = []
        residues = set()

        for protein_id in protein_ids:
            protein_residue = protein_mol[protein_id]
            metadata = fp.metadata(ligand_residue, protein_residue)

            for interaction_name, occurrences in (metadata or {}).items():
                residues.add(protein_id)
                for occurrence in occurrences:
                    records.append(
                        {
                            "ligand_residue": str(ligand_residue.resid),
                            "protein_residue": protein_id,
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

    summary_data = []
    for lid in ligands:
        records = all_results[lid]
        residues = ligand_residues[lid]
        interaction_types = Counter(
            record["interaction"] for record in records
        )

        affinity = float(
            ranking.loc[
                ranking["ligand_id"] == lid, "affinity"
            ].iloc[0]
        )

        summary_data.append(
            {
                "ligand_id": lid,
                "affinity": affinity,
                "n_interactions": len(records),
                "n_residues": len(residues),
                "interaction_types": "; ".join(
                    f"{name}:{count}"
                    for name, count in interaction_types.most_common()
                ),
            }
        )

    summary = (
        pd.DataFrame(summary_data)
        .sort_values("affinity")
        .reset_index(drop=True)
    )

    if len(summary) != TOP_N or set(summary["ligand_id"]) != set(ligands):
        raise RuntimeError("ProLIF summary does not contain exactly the current Top 5")

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, index=False)

    print(f"\n   Saved: {output}")


if __name__ == "__main__":
    main()
