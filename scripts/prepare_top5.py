import os
import subprocess
import tempfile

import pandas as pd
from rdkit import Chem

DOCKING_RESULTS = "docking/results/standardized_affinities.csv"
POSE_DIR = "docking/results/standardized"
OUTPUT_DIR = "docking/results/fixed"

df_all = pd.read_csv(DOCKING_RESULTS)
df_all["affinity"] = pd.to_numeric(df_all["affinity"], errors="coerce")
df_all = df_all[df_all["status"].isin(["ok", "existing"])].dropna(subset=["affinity"])
df_all = df_all.sort_values("affinity").reset_index(drop=True)
top_5 = df_all.head(5)["ligand_id"].tolist()

print(f"Standardized Top 5: {top_5}\n")

if len(top_5) != 5:
    raise RuntimeError(f"Expected 5 Top hits, found {len(top_5)}")

os.makedirs(OUTPUT_DIR, exist_ok=True)


def remove_all_hydrogens(mol):
    """Remove every H atom, including isolated/malformed H atoms."""
    editable = Chem.RWMol(mol)
    h_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.GetSymbol() == "H"
    ]
    for idx in sorted(h_indices, reverse=True):
        editable.RemoveAtom(idx)
    return editable.GetMol()


def element_counts(mol):
    heavy = remove_all_hydrogens(mol)
    return sorted(atom.GetSymbol() for atom in heavy.GetAtoms())


def convert_pdbqt_to_sdf(pdbqt_path, output_sdf):
    """
    Convert the standardized PDBQT pose directly to SDF.

    This intentionally does not force atom correspondence with the original
    SDF template: the docking/preparation pipeline may reorder atoms. The
    standardized PDBQT is therefore the authoritative pose representation.
    """
    result = subprocess.run(
        [
            "obabel",
            "-ipdbqt",
            pdbqt_path,
            "-osdf",
            "-O",
            output_sdf,
        ],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0 or not os.path.exists(output_sdf):
        detail = (result.stderr or result.stdout).strip()
        return False, f"Open Babel failed: {detail[:160]}"

    return True, "PDBQT -> SDF"


def validate_pose(ligand_id, pose_sdf):
    """Validate that the converted pose is chemically consistent with its input ligand."""
    template_path = f"ligands/{ligand_id}.sdf"

    template = Chem.MolFromMolFile(
        template_path,
        removeHs=False,
        sanitize=False,
    )
    pose = Chem.MolFromMolFile(
        pose_sdf,
        removeHs=False,
        sanitize=False,
    )

    if template is None:
        return False, f"Cannot read template {template_path}"
    if pose is None:
        return False, "Cannot read converted pose"

    template_symbols = element_counts(template)
    pose_symbols = element_counts(pose)

    if template_symbols != pose_symbols:
        return (
            False,
            "Element composition mismatch "
            f"(template={template_symbols}, pose={pose_symbols})",
        )

    try:
        Chem.SanitizeMol(pose)
    except Exception as exc:
        return False, f"RDKit sanitization failed: {exc}"

    return True, "validated"


success, failed = [], []

print("Preparing standardized Top 5 for ProLIF...\n")

for lid in top_5:
    pdbqt = f"{POSE_DIR}/{lid}_out.pdbqt"
    output = f"{OUTPUT_DIR}/{lid}_fixed.sd"

    if not os.path.exists(pdbqt):
        print(f"  FAIL {lid:5s} — standardized pose not found")
        failed.append(lid)
        continue

    ok, method = convert_pdbqt_to_sdf(pdbqt, output)

    if not ok:
        print(f"  FAIL {lid:5s} — {method}")
        failed.append(lid)
        continue

    ok, validation = validate_pose(lid, output)

    if ok:
        print(f"  OK   {lid:5s} via {method}; {validation}")
        success.append(lid)
    else:
        print(f"  FAIL {lid:5s} — {validation}")
        failed.append(lid)

print(f"\n{'=' * 60}")
print(f"  Success: {len(success)}/{len(top_5)}")
print(f"  Failed:  {len(failed)}/{len(top_5)}")
print(f"{'=' * 60}")

if failed:
    print(f"\nFailed ligands: {failed}")
