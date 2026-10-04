"""
merge_final_results.py — Merge Docking, ADMET, and ProLIF results into a final summary.

This script combines results from different stages of the drug discovery
pipeline into a single ranked table of candidate molecules.

Usage
-----
    python scripts/merge_final_results.py

Output
------
    results/top_hits_summary.csv
"""

import sys
from pathlib import Path

import pandas as pd


def load_docking_results(path: str = "docking/results/affinities_recovered_111.csv") -> pd.DataFrame:
    """Load and clean docking affinity results.

    Parameters
    ----------
    path : str
        Path to the affinities CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ``ligand_id``, ``affinity``, ``rank``.
    """
    if not Path(path).exists():
        print(f"  [WARN] Docking results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)
    print(f"  [OK] Docking results: {len(df)} ligands loaded")
    return df


def load_prolif_results(path: str = "docking/results/prolif_summary.csv") -> pd.DataFrame:
    """Load ProLIF summary results.

    Parameters
    ----------
    path : str
        Path to ``prolif_summary.csv``.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ``ligand_id``, ``n_interactions``.
    """
    results_path = Path(path)

    if not results_path.exists():
        print(f"  [WARN] ProLIF results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(results_path)

    required = {"ligand_id", "n_interactions"}
    missing = required - set(df.columns)

    if missing:
        print(f"  [WARN] ProLIF summary missing columns: {sorted(missing)}")
        return pd.DataFrame()

    result = df[["ligand_id", "n_interactions"]].copy()
    print(f"  [OK] ProLIF results: {len(result)} ligands loaded")
    return result


def merge_results(df_docking: pd.DataFrame, df_prolif: pd.DataFrame) -> pd.DataFrame:
    """Merge docking and ProLIF results.

    Parameters
    ----------
    df_docking : pd.DataFrame
        Docking affinities.
    df_prolif : pd.DataFrame
        ProLIF interaction counts.

    Returns
    -------
    pd.DataFrame
        Merged and sorted DataFrame.
    """
    if df_docking.empty:
        return pd.DataFrame()

    df = df_docking.copy()

    if not df_prolif.empty:
        df = df.merge(df_prolif, on="ligand_id", how="left")
        df["n_interactions"] = pd.to_numeric(
            df["n_interactions"],
            errors="coerce",
        )

    df = df.sort_values("affinity")
    return df


def main():
    """Main entry point — merge all results and save summary."""
    print("=" * 60)
    print("  Merging Final Results")
    print("=" * 60)

    # Load results
    print("\n1. Loading results...")
    df_docking = load_docking_results()
    df_prolif = load_prolif_results()

    # Merge
    print("\n2. Merging...")
    df_final = merge_results(df_docking, df_prolif)

    if df_final.empty:
        print("\n  [ERROR] No results to merge. Run docking first.")
        sys.exit(1)

    # Save
    output_dir = Path("results")
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "top_hits_summary.csv"
    df_final.to_csv(output_path, index=False)

    # Print summary
    print("\n3. Summary (Top 10):")
    print("-" * 60)
    print(df_final.head(10).to_string(index=False))
    print(f"\n  Total ligands: {len(df_final)}")
    print(f"  Saved: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
