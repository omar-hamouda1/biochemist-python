from pathlib import Path

import pandas as pd


OLD_RESULTS = Path("docking/results/affinities_full.csv")
RECOVERED_RESULTS = Path("docking/results/recovered_affinities.csv")
SDF_DIR = Path("ligands")

OUTPUT = Path("docking/results/affinities_recovered_111.csv")


def main():
    old = pd.read_csv(OLD_RESULTS)
    recovered = pd.read_csv(RECOVERED_RESULTS)

    old = old[["ligand_id", "affinity"]].copy()
    recovered = recovered[["ligand_id", "affinity"]].copy()

    old_ids = set(old["ligand_id"])
    recovered_ids = set(recovered["ligand_id"])

    duplicates = sorted(old_ids & recovered_ids)

    if duplicates:
        raise ValueError(
            f"Duplicate ligand IDs found: {duplicates}"
        )

    sdf_ids = {path.stem for path in SDF_DIR.glob("*.sdf")}

    unknown_ids = sorted(
        (old_ids | recovered_ids) - sdf_ids
    )

    if unknown_ids:
        raise ValueError(
            f"Results contain unknown ligands: {unknown_ids}"
        )

    merged = pd.concat(
        [old, recovered],
        ignore_index=True,
    )

    merged["affinity"] = pd.to_numeric(
        merged["affinity"],
        errors="raise",
    )

    merged = (
        merged
        .sort_values("affinity", ascending=True)
        .reset_index(drop=True)
    )

    merged.insert(
        0,
        "rank",
        range(1, len(merged) + 1),
    )

    expected = len(old) + len(recovered)

    if len(merged) != expected:
        raise ValueError(
            f"Expected {expected} results, got {len(merged)}"
        )

    merged.to_csv(OUTPUT, index=False)

    print("=" * 70)
    print("DOCKING RECOVERY MERGE")
    print("=" * 70)
    print(f"Original results : {len(old)}")
    print(f"Recovered results: {len(recovered)}")
    print(f"Final unique     : {len(merged)}")
    print(f"Missing ligands  : {len(sdf_ids - set(merged['ligand_id']))}")
    print(f"Output           : {OUTPUT}")
    print("=" * 70)


if __name__ == "__main__":
    main()
