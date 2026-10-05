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


def convert_pose_to_pdb(ligand_id, pdbqt_pose):
    """Convert PDBQT to PDB without adding hydrogens."""
    temp_pdb = os.path.join(
        tempfile.gettempdir(), f"{ligand_id}_standardized_raw.pdb"
    )
    result = subprocess.run(
        ["obabel", "-ipdbqt", pdbqt_pose, "-opdb", "-O", temp_pdb],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not os.path.exists(temp_pdb):
        return None
    return temp_pdb


def remove_all_hydrogens(mol):
    """Remove every H atom by index, including malformed/isolated H atoms."""
    editable = Chem.RWMol(mol)
    h_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.GetSymbol() == "H"
    ]
    for idx in sorted(h_indices, reverse=True):
        editable.RemoveAtom(idx)
    return editable.GetMol()


def method_template_coordinates(temp_pdb, template_sdf, output_sdf):
    """
    Preserve bond orders/connectivity from the original ligand template and
    transfer docked heavy-atom coordinates from the standardized pose.
    """
    try:
        pose = Chem.MolFromPDBFile(
            temp_pdb, removeHs=True, sanitize=False
        )
        template = Chem.MolFromMolFile(
            template_sdf, removeHs=False, sanitize=False
        )

        if pose is None:
            return False, "Cannot read docked PDB"
        if template is None:
            return False, "Cannot read template"

        template_noH = remove_all_hydrogens(template)

        if pose.GetNumAtoms() != template_noH.GetNumAtoms():
            return (
                False,
                f"Heavy-atom mismatch ({pose.GetNumAtoms()} vs "
                f"{template_noH.GetNumAtoms()})",
            )

        pose_symbols = [a.GetSymbol() for a in pose.GetAtoms()]
        template_symbols = [a.GetSymbol() for a in template_noH.GetAtoms()]
        if pose_symbols != template_symbols:
            return False, "Heavy-atom order/type mismatch"

        pose_conf = pose.GetConformer()
        final_mol = Chem.Mol(template_noH)
        final_mol.RemoveAllConformers()

        conf = Chem.Conformer(template_noH.GetNumAtoms())
        for idx in range(template_noH.GetNumAtoms()):
            conf.SetAtomPosition(idx, pose_conf.GetAtomPosition(idx))
        final_mol.AddConformer(conf, assignId=True)

        final_mol = Chem.AddHs(final_mol, addCoords=True)
        Chem.SanitizeMol(final_mol)

        writer = Chem.SDWriter(output_sdf)
        writer.write(final_mol)
        writer.close()
        return True, "template + docked coordinates"
    except Exception as exc:
        return False, f"template-coordinate error: {exc}"


def method_auto(temp_pdb, output_sdf):
    """Last-resort bond inference from the docked pose."""
    try:
        mol = Chem.MolFromPDBFile(
            temp_pdb, removeHs=True, sanitize=False
        )
        if mol is None:
            return False, "Cannot read docked PDB"

        from rdkit.Chem import rdDetermineBonds

        for charge in [0, 1, -1, 2, -2]:
            try:
                candidate = Chem.Mol(mol)
                rdDetermineBonds.DetermineBonds(candidate, charge=charge)
                candidate = Chem.AddHs(candidate, addCoords=True)
                Chem.SanitizeMol(candidate)

                writer = Chem.SDWriter(output_sdf)
                writer.write(candidate)
                writer.close()
                return True, f"auto (q={charge})"
            except Exception:
                continue

        return False, "all charges failed"
    except Exception as exc:
        return False, f"auto error: {exc}"


success, failed = [], []
print("Preparing standardized Top 5 for ProLIF...\n")

for lid in top_5:
    pdbqt = f"{POSE_DIR}/{lid}_out.pdbqt"
    output = f"{OUTPUT_DIR}/{lid}_fixed.sd"

    if not os.path.exists(pdbqt):
        print(f"  FAIL {lid}: standardized pose not found")
        failed.append(lid)
        continue

    temp_pdb = convert_pose_to_pdb(lid, pdbqt)
    if temp_pdb is None:
        print(f"  FAIL {lid}: Open Babel conversion failed")
        failed.append(lid)
        continue

    template = "pdb/ligand.sd" if lid == "13U" else f"ligands/{lid}.sdf"

    ok, method = method_template_coordinates(
        temp_pdb, template, output
    )

    if not ok:
        ok, method = method_auto(temp_pdb, output)

    if ok:
        print(f"  OK   {lid:5s} via {method}")
        success.append(lid)
    else:
        print(f"  FAIL {lid:5s} — {method}")
        failed.append(lid)

print(f"\n{'=' * 60}")
print(f"  Success: {len(success)}/{len(top_5)}")
print(f"  Failed:  {len(failed)}/{len(top_5)}")
print(f"{'=' * 60}")

if failed:
    print(f"\nFailed ligands: {failed}")
