"""Merge standardized docking, ProLIF, and ADMET results."""

import sys
from pathlib import Path

import pandas as pd


DOCKING_PATH = "docking/results/standardized_affinities.csv"
PROLIF_PATH = "docking/results/prolif_summary.csv"
ADMET_PATH = "docking/results/admet/admet_top5.csv"
OUTPUT_PATH = "results/top_hits_summary.csv"


def load_docking_results(path=DOCKING_PATH):
    if not Path(path).exists():
        print(f"  [ERROR] Docking results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    required = {"ligand_id", "affinity", "status"}
    missing = required - set(df.columns)
    if missing:
        print(
            "  [ERROR] Docking results missing columns: "
            f"{sorted(missing)}"
        )
        return pd.DataFrame()

    df = df[df["status"].isin(["ok", "existing"])].copy()

    df["affinity"] = pd.to_numeric(
        df["affinity"],
        errors="coerce",
    )

    df = (
        df.dropna(subset=["affinity"])
        .sort_values("affinity")
        .reset_index(drop=True)
    )

    if df["ligand_id"].duplicated().any():
        print("  [ERROR] Duplicate ligand IDs in docking results")
        return pd.DataFrame()

    df["rank"] = df.index + 1

    print(
        f"  [OK] Standardized docking: "
        f"{len(df)} valid ligands loaded"
    )

    return df


def load_prolif_results(path=PROLIF_PATH):
    if not Path(path).exists():
        print(f"  [INFO] ProLIF results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    required = {"ligand_id", "n_interactions"}
    missing = required - set(df.columns)

    if missing:
        print(
            "  [WARN] ProLIF summary missing columns: "
            f"{sorted(missing)}"
        )
        return pd.DataFrame()

    if df["ligand_id"].duplicated().any():
        print("  [WARN] Duplicate ligand IDs in ProLIF results")
        return pd.DataFrame()

    # Keep ProLIF-specific columns only.
    keep = ["ligand_id"]

    for column in [
        "n_interactions",
        "n_residues",
        "interaction_types",
    ]:
        if column in df.columns:
            keep.append(column)

    result = df[keep].copy()

    result["n_interactions"] = pd.to_numeric(
        result["n_interactions"],
        errors="coerce",
    )

    if "n_residues" in result.columns:
        result["n_residues"] = pd.to_numeric(
            result["n_residues"],
            errors="coerce",
        )

    print(
        f"  [OK] ProLIF results: "
        f"{len(result)} ligands loaded"
    )

    return result


def load_admet_results(path=ADMET_PATH):
    if not Path(path).exists():
        print(f"  [INFO] ADMET results not found: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path)

    required = {
        "ligand_id",
        "Lipinski_Verdict",
        "Veber_Verdict",
        "PAINS_Verdict",
    }

    missing = required - set(df.columns)

    if missing:
        print(
            "  [WARN] ADMET summary missing columns: "
            f"{sorted(missing)}"
        )
        return pd.DataFrame()

    if df["ligand_id"].duplicated().any():
        print("  [WARN] Duplicate ligand IDs in ADMET results")
        return pd.DataFrame()

    # Docking affinity is the authoritative affinity.
    # Remove duplicate affinity from ADMET before merging.
    keep = [
        column
        for column in df.columns
        if column not in {"ligand_id", "affinity"}
    ]

    result = df[["ligand_id", *keep]].copy()

    print(
        f"  [OK] ADMET results: "
        f"{len(result)} ligands loaded"
    )

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

    print("\n3. Loading ADMET...")
    admet = load_admet_results()

    print("\n4. Merging...")

    final = docking.copy()

    if not prolif.empty:
        final = final.merge(
            prolif,
            on="ligand_id",
            how="left",
            validate="one_to_one",
        )

    if not admet.empty:
        final = final.merge(
            admet,
            on="ligand_id",
            how="left",
            validate="one_to_one",
        )

    # Final schema must contain exactly one authoritative affinity column.
    if "affinity" not in final.columns:
        print("  [ERROR] Final result lost authoritative docking affinity")
        sys.exit(1)

    forbidden = [
        column
        for column in final.columns
        if column in {"affinity_x", "affinity_y"}
    ]

    if forbidden:
        print(
            "  [ERROR] Unexpected duplicate affinity columns: "
            f"{forbidden}"
        )
        sys.exit(1)

    # Top-5 consistency check across downstream analyses.
    docking_top5 = docking.head(5)["ligand_id"].tolist()

    if not prolif.empty:
        prolif_ids = set(prolif["ligand_id"])
        missing = [ligand for ligand in docking_top5 if ligand not in prolif_ids]
        if missing:
            print(
                "  [ERROR] Top docking hits missing from ProLIF: "
                f"{missing}"
            )
            sys.exit(1)

    if not admet.empty:
        admet_ids = set(admet["ligand_id"])
        missing = [ligand for ligand in docking_top5 if ligand not in admet_ids]
        if missing:
            print(
                "  [ERROR] Top docking hits missing from ADMET: "
                f"{missing}"
            )
            sys.exit(1)

    Path(OUTPUT_PATH).parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    final.to_csv(
        OUTPUT_PATH,
        index=False,
    )

    display_columns = [
        "rank",
        "ligand_id",
        "affinity",
        "n_interactions",
        "n_residues",
        "Lipinski_Verdict",
        "Veber_Verdict",
        "PAINS_Verdict",
    ]

    available = [
        column
        for column in display_columns
        if column in final.columns
    ]

    print("\n5. Standardized Top 10:")
    print("-" * 60)

    print(
        final.head(10)[available].to_string(index=False)
    )

    print(f"\n  Total ligands: {len(final)}")
    print(f"  Top 5: {docking_top5}")
    print(f"  Saved: {OUTPUT_PATH}")
    print("=" * 60)


if __name__ == "__main__":
    main()
