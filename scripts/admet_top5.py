import os

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, FilterCatalog, FilterCatalogParams, Lipinski, rdMolDescriptors

DOCKING_RESULTS = "docking/results/standardized_affinities.csv"
FIXED_DIR = "docking/results/fixed"
OUTPUT_DIR = "docking/results/admet"
TOP_N = 5

os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_top_hits():
    ranking = pd.read_csv(DOCKING_RESULTS)
    ranking["affinity"] = pd.to_numeric(ranking["affinity"], errors="coerce")
    ranking = ranking[
        ranking["status"].isin(["ok", "existing"])
    ].dropna(subset=["affinity"])
    return ranking.sort_values("affinity").head(TOP_N).reset_index(drop=True)


def load_ligand(ligand_id):
    path = os.path.join(FIXED_DIR, f"{ligand_id}_fixed.sd")
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=True)
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
    return FilterCatalog.FilterCatalog(params)


ranking = load_top_hits()

print("=" * 72)
print("  ADMET / Drug-likeness Analysis: Standardized Top 5")
print("=" * 72)
print("\nTop 5:")
for rank, row in ranking.iterrows():
    print(f"  {rank + 1}. {row['ligand_id']:5s}  {row['affinity']:.6f}")

catalog = pains_catalog()
rows = []

for _, row in ranking.iterrows():
    lid = row["ligand_id"]
    mol, path = load_ligand(lid)
    props = descriptors(mol)
    lip_violations, lip_verdict = lipinski(props)
    veb_violations, veb_verdict = veber(props)

    matches = catalog.GetMatches(mol)
    alerts = [match.GetDescription() for match in matches]

    rows.append(
        {
            "rank": int(len(rows) + 1),
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
            "source_pose": path,
        }
    )

summary = pd.DataFrame(rows)

summary.to_csv(
    os.path.join(OUTPUT_DIR, "admet_summary.csv"),
    index=False,
)

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
].to_csv(
    os.path.join(OUTPUT_DIR, "admet_top5.csv"),
    index=False,
)

print("\n" + "=" * 72)
print("  Final ADMET / Drug-likeness Summary")
print("=" * 72)
print(
    summary[
        [
            "rank",
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
    ].to_string(index=False)
)

print(f"\nSaved: {OUTPUT_DIR}/admet_summary.csv")
print(f"Saved: {OUTPUT_DIR}/admet_top5.csv")
print("\nNote: these are RDKit descriptors and rule-based filters, not full ADME/Toxicity predictions.")
