import os
import subprocess
import tempfile

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

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
    temp_pdb = os.path.join(tempfile.gettempdir(), f"{ligand_id}_standardized_raw.pdb")
    result = subprocess.run(["obabel", "-ipdbqt", pdbqt_pose, "-opdb", "-O", temp_pdb, "-h"], capture_output=True, text=True)
    if result.returncode != 0:
        return None
    return temp_pdb if os.path.exists(temp_pdb) else None


def method_template(temp_pdb, template_sdf, output_sdf):
    try:
        pose = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
        if pose is None:
            return False, "Cannot read PDB"
        pose_noH = Chem.RemoveHs(pose, sanitize=False)
        template = Chem.MolFromMolFile(template_sdf, removeHs=False)
        if template is None:
            return False, "Cannot read template"
        template_noH = Chem.RemoveHs(template)
        remaining_h = [a.GetIdx() for a in template_noH.GetAtoms() if a.GetSymbol() == "H"]
        if remaining_h:
            editable = Chem.RWMol(template_noH)
            for idx in sorted(remaining_h, reverse=True):
                editable.RemoveAtom(idx)
            template_noH = editable.GetMol()
        if pose_noH.GetNumAtoms() != template_noH.GetNumAtoms():
            return False, "Atom mismatch"
        final_mol = AllChem.AssignBondOrdersFromTemplate(template_noH, pose_noH)
        final_mol = Chem.AddHs(final_mol)
        Chem.SanitizeMol(final_mol)
        writer = Chem.SDWriter(output_sdf)
        writer.write(final_mol)
        writer.close()
        return True, "template"
    except Exception as exc:
        return False, str(exc)[:80]


def method_auto(temp_pdb, output_sdf):
    for charge in [0, 1, -1, 2, -2]:
        try:
            mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
            if mol is None:
                continue
            mol_noH = Chem.RemoveHs(mol)
            rdDetermineBonds.DetermineBonds(mol_noH, charge=charge)
            mol_final = Chem.AddHs(mol_noH)
            Chem.SanitizeMol(mol_final)
            writer = Chem.SDWriter(output_sdf)
            writer.write(mol_final)
            writer.close()
            return True, f"auto (q={charge})"
        except Exception:
            continue
    return False, "all charges failed"

success, failed = [], []
print("Preparing standardized Top 5 for ProLIF...\n")

for lid in top_5:
    pdbqt = f"{POSE_DIR}/{lid}_out.pdbqt"
    output_sdf = f"{OUTPUT_DIR}/{lid}_fixed.sd"
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
    ok, method = method_template(temp_pdb, template, output_sdf)
    if not ok:
        ok, method = method_auto(temp_pdb, output_sdf)
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
