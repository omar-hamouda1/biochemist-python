"""Merge current standardized docking, ProLIF, and ADMET results."""

from pathlib import Path
import sys

import pandas as pd

from src.docking_config import load_docking_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TOP_N = 5
PROLIF_PATH = PROJECT_ROOT / "docking" / "results" / "prolif_summary.csv"
ADMET_PATH = PROJECT_ROOT / "docking" / "results" / "admet" / "admet_top5.csv"
OUTPUT_PATH = PROJECT_ROOT / "results" / "top_hits_summary.csv"


def load_docking_results(path: Path, expected_ligands: int):
    if not path.exists():
        print(f"  [ERROR] Docking results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    required = {"ligand_id", "affinity", "status"}
    missing = required - set(df.columns)
    if missing:
        print(f"  [ERROR] Docking results missing columns: {sorted(missing)}")
        return pd.DataFrame()

    df = df[df["status"].isin(["ok", "existing"])].copy()
    df["affinity"] = pd.to_numeric(df["affinity"], errors="coerce")
    df = (
        df.dropna(subset=["affinity"])
        .sort_values("affinity")
        .reset_index(drop=True)
    )

    if len(df) != expected_ligands:
        print(
            f"  [ERROR] Expected {expected_ligands} valid docking rows, "
            f"found {len(df)}"
        )
        return pd.DataFrame()

    if df["ligand_id"].duplicated().any():
        print("  [ERROR] Duplicate ligand IDs in docking results")
        return pd.DataFrame()

    df["rank"] = df.index + 1
    print(f"  [OK] Standardized docking: {len(df)} valid ligands loaded")
    return df


def load_downstream_results(path: Path, required_columns, label: str):
    if not path.exists():
        print(f"  [ERROR] {label} results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    missing = set(required_columns) - set(df.columns)
    if missing:
        print(f"  [ERROR] {label} missing columns: {sorted(missing)}")
        return pd.DataFrame()

    if df["ligand_id"].duplicated().any():
        print(f"  [ERROR] Duplicate ligand IDs in {label} results")
        return pd.DataFrame()

    return df


def main():
    config = load_docking_config()

    print("=" * 60)
    print("  Merging Current Standardized Results")
    print("=" * 60)

    print("\n1. Loading standardized docking...")
    docking = load_docking_results(config.report, config.expected_ligands)
    if docking.empty:
        sys.exit(1)

    docking_top5 = docking.head(TOP_N)["ligand_id"].astype(str).tolist()

    print("\n2. Loading current ProLIF...")
    prolif = load_downstream_results(
        PROLIF_PATH,
        {"ligand_id", "n_interactions"},
        "ProLIF",
    )
    if prolif.empty:
        sys.exit(1)

    print("\n3. Loading current rule-based ADMET...")
    admet = load_downstream_results(
        ADMET_PATH,
        {"ligand_id", "Lipinski_Verdict", "Veber_Verdict", "PAINS_Verdict"},
        "ADMET",
    )
    if admet.empty:
        sys.exit(1)

    for label, df in (("ProLIF", prolif), ("ADMET", admet)):
        ids = set(df["ligand_id"].astype(str))
        missing = [ligand for ligand in docking_top5 if ligand not in ids]
        if missing:
            print(f"  [ERROR] Current Top 5 missing from {label}: {missing}")
            sys.exit(1)

    final = docking.copy()

    prolif_keep = [
        column
        for column in [
            "ligand_id",
            "n_interactions",
            "n_residues",
            "interaction_types",
        ]
        if column in prolif.columns
    ]
    final = final.merge(
        prolif[prolif_keep],
        on="ligand_id",
        how="left",
        validate="one_to_one",
    )

    admet_keep = [
        column for column in admet.columns
        if column not in {"ligand_id", "affinity"}
    ]
    final = final.merge(
        admet[["ligand_id", *admet_keep]],
        on="ligand_id",
        how="left",
        validate="one_to_one",
    )

    if "affinity" not in final.columns:
        print("  [ERROR] Final result lost authoritative docking affinity")
        sys.exit(1)

    forbidden = [
        column for column in final.columns
        if column in {"affinity_x", "affinity_y"}
    ]
    if forbidden:
        print(f"  [ERROR] Unexpected duplicate affinity columns: {forbidden}")
        sys.exit(1)

    Path(OUTPUT_PATH).parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUTPUT_PATH, index=False)

    print("\n4. Standardized Top 10:")
    print(
        final.head(10).to_string(index=False)
    )
    print(f"\n  Total ligands: {len(final)}")
    print(f"  Top 5: {docking_top5}")
    print(f"  Saved: {OUTPUT_PATH}")
    print("=" * 60)


if __name__ == "__main__":
    main()
