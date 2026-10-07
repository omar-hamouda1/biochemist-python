"""Calculate rule-based molecular-property filters for the current Top 5."""

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

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
    ranking = ranking.sort_values(["affinity", "ligand_id"], kind="mergesort").reset_index(drop=True)

    top = ranking.head(TOP_N)
    if len(top) != TOP_N or top["ligand_id"].duplicated().any():
        raise RuntimeError("Current standardized Top 5 is incomplete or duplicated")

    return top


def load_ligand(ligand_id, fixed_dir):
    path = fixed_dir / f"{ligand_id}_fixed.sdf"
    if not path.exists():
        raise FileNotFoundError(path)

    mol = Chem.MolFromMolFile(
        str(path), removeHs=False, sanitize=True
    )
    if mol is None:
        raise ValueError(f"RDKit could not parse {path}")

    return mol, path


def descriptors(mol):
    return {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "RotB": Lipinski.NumRotatableBonds(mol),
        "Rings": Lipinski.RingCount(mol),
        "Aromatic_Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Heavy_Atoms": mol.GetNumHeavyAtoms(),
        "Fraction_Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Mol_Refract": round(Descriptors.MolMR(mol), 2),
        "Formal_Charge": Chem.GetFormalCharge(mol),
    }


def lipinski(props):
    checks = {
        "MW": props["MW"] <= 500,
        "LogP": props["LogP"] <= 5,
        "HBD": props["HBD"] <= 5,
        "HBA": props["HBA"] <= 10,
    }
    violations = sum(not passed for passed in checks.values())
    verdict = (
        "Excellent"
        if violations == 0
        else "Acceptable"
        if violations == 1
        else "Poor"
    )
    return violations, verdict


def veber(props):
    violations = int(props["RotB"] > 10) + int(props["TPSA"] > 140)
    return violations, ("Pass" if violations == 0 else "Fail")


def pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def main():
    config = load_docking_config()
    ranking = load_top_hits(config.report)
    fixed_dir = PROJECT_ROOT / "docking" / "results" / "fixed"
    output_dir = PROJECT_ROOT / "docking" / "results" / "admet"
    output_dir.mkdir(parents=True, exist_ok=True)

    catalog = pains_catalog()
    rows = []

    print("=" * 72)
    print("  ADMET / Drug-likeness Analysis: Standardized Top 5")
    print("=" * 72)

    for rank, (_, row) in enumerate(ranking.iterrows(), start=1):
        lid = str(row["ligand_id"])
        mol, path = load_ligand(lid, fixed_dir)
        props = descriptors(mol)
        lip_violations, lip_verdict = lipinski(props)
        veb_violations, veb_verdict = veber(props)

        matches = catalog.GetMatches(mol)
        alerts = [match.GetDescription() for match in matches]

        rows.append(
            {
                "rank": rank,
                "ligand_id": lid,
                "affinity": float(row["affinity"]),
                **props,
                "Lipinski_Violations": lip_violations,
                "Lipinski_Verdict": lip_verdict,
                "Veber_Violations": veb_violations,
                "Veber_Verdict": veb_verdict,
                "PAINS_Alerts": len(alerts),
                "PAINS_Verdict": "Clean" if not alerts else "Flagged",
                "PAINS_Details": "; ".join(alerts),
                "source_pose": str(path.relative_to(PROJECT_ROOT)),
            }
        )

    summary = pd.DataFrame(rows)

    if len(summary) != TOP_N or set(summary["ligand_id"]) != set(ranking["ligand_id"]):
        raise RuntimeError("ADMET summary does not contain exactly the current Top 5")

    summary.to_csv(output_dir / "admet_summary.csv", index=False)

    summary[
        [
            "ligand_id",
            "affinity",
            "MW",
            "LogP",
            "HBD",
            "HBA",
            "TPSA",
            "RotB",
            "Lipinski_Violations",
            "Lipinski_Verdict",
            "Veber_Violations",
            "Veber_Verdict",
            "PAINS_Alerts",
            "PAINS_Verdict",
        ]
    ].to_csv(output_dir / "admet_top5.csv", index=False)

    print(summary.to_string(index=False))
    print(f"\nSaved: {output_dir / 'admet_summary.csv'}")
    print(f"Saved: {output_dir / 'admet_top5.csv'}")
    print(
        "\nNote: these are RDKit descriptors and rule-based filters, "
        "not full ADME/Toxicity predictions."
    )


if __name__ == "__main__":
    main()
