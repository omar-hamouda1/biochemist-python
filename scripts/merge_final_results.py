"""Merge authoritative standardized docking with available ProLIF results."""

import sys
from pathlib import Path
import pandas as pd

DOCKING_PATH = "docking/results/standardized_affinities.csv"
PROLIF_PATH = "docking/results/prolif_summary.csv"
OUTPUT_PATH = "results/top_hits_summary.csv"


def load_docking_results(path=DOCKING_PATH):
    if not Path(path).exists():
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
    df = df.dropna(subset=["affinity"]).sort_values("affinity").reset_index(drop=True)
    df["rank"] = df.index + 1
    print(f"  [OK] Standardized docking: {len(df)} valid ligands loaded")
    return df


def load_prolif_results(path=PROLIF_PATH):
    if not Path(path).exists():
        print(f"  [INFO] ProLIF results not found: {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    required = {"ligand_id", "n_interactions"}
    missing = required - set(df.columns)
    if missing:
        print(f"  [WARN] ProLIF summary missing columns: {sorted(missing)}")
        return pd.DataFrame()
    result = df[["ligand_id", "n_interactions"]].copy()
    result["n_interactions"] = pd.to_numeric(result["n_interactions"], errors="coerce")
    print(f"  [OK] ProLIF results: {len(result)} ligands loaded")
    return result


def main():
    print("=" * 60)
    print("  Merging Final Results")
    print("=" * 60)
    print("\n1. Loading standardized docking...")
    docking = load_docking_results()
    if docking.empty:
        sys.exit(1)
    print("\n2. Loading ProLIF...")
    prolif = load_prolif_results()
    print("\n3. Merging...")
    final = docking.merge(prolif, on="ligand_id", how="left") if not prolif.empty else docking
    Path(OUTPUT_PATH).parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUTPUT_PATH, index=False)
    print("\n4. Standardized Top 10:")
    print("-" * 60)
    print(final.head(10)[["rank", "ligand_id", "affinity", "n_interactions"]].to_string(index=False))
    print(f"\n  Total ligands: {len(final)}")
    print(f"  Saved: {OUTPUT_PATH}")
    print("=" * 60)


if __name__ == "__main__":
    main()
