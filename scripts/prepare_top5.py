import os
import subprocess
import tempfile
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

# ─── Load docking results ───
df_all = pd.read_csv("docking/results/affinities_recovered_111.csv")
top_5 = df_all.head(5)["ligand_id"].tolist()
print(f"Top 5 ligands: {top_5}\n")

output_dir = "docking/results/fixed"
os.makedirs(output_dir, exist_ok=True)


def convert_pose_to_pdb(ligand_id, pdbqt_pose):
    """PDBQT → PDB using obabel."""
    temp_pdb = os.path.join(tempfile.gettempdir(), f"{ligand_id}_raw.pdb")
    subprocess.run(
        ["obabel", "-ipdbqt", pdbqt_pose, "-opdb", "-O", temp_pdb, "-h"], capture_output=True
    )
    return temp_pdb if os.path.exists(temp_pdb) else None


def method_template(ligand_id, temp_pdb, template_sdf, output_sdf):
    """Template-based method."""
    try:
        pose = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
        if pose is None:
            return False, "Cannot read PDB"

        pose_noH = Chem.RemoveHs(pose, sanitize=False)

        template = Chem.MolFromMolFile(template_sdf, removeHs=False)
        if template is None:
            return False, "Cannot read template"

        template_noH = Chem.RemoveHs(template)

        remaining_h = [
            atom.GetIdx()
            for atom in template_noH.GetAtoms()
            if atom.GetSymbol() == "H"
        ]

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
    except Exception as e:
        return False, str(e)[:50]


def method_auto(ligand_id, temp_pdb, output_sdf):
    """Auto-detect method (try charges 0, ±1, ±2)."""
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


# ─── Process each ligand (hybrid) ───
print("Processing top 5 with hybrid method...\n")
success, failed = [], []

for lid in top_5:
    pdbqt = f"docking/results/{lid}_out.pdbqt"

    if not os.path.exists(pdbqt):
        print(f"  FAIL {lid}: pose not found")
        failed.append(lid)
        continue

    # Convert PDBQT → PDB once
    temp_pdb = convert_pose_to_pdb(lid, pdbqt)
    if temp_pdb is None:
        print(f"  FAIL {lid}: obabel failed")
        failed.append(lid)
        continue

    output_sdf = f"{output_dir}/{lid}_fixed.sd"

    # Determine template
    if lid == "13U":
        template = "pdb/ligand.sd"
    else:
        template = f"ligands/{lid}.sdf"

    # ─── Try 1: Template method ───
    ok, method = method_template(lid, temp_pdb, template, output_sdf)

    # ─── Try 2: Auto method (if template failed) ───
    if not ok:
        ok, method = method_auto(lid, temp_pdb, output_sdf)

    if ok:
        print(f"  OK   {lid:5s} via {method}")
        success.append(lid)
    else:
        print(f"  FAIL {lid:5s} — both methods failed")
        failed.append(lid)

print(f"\n{'='*60}")
print(f"  Success: {len(success)}/{len(top_5)}")
print(f"  Failed:  {len(failed)}/{len(top_5)}")
print(f"{'='*60}")

if failed:
    print(f"\nFailed ligands: {failed}")
    print(f"Note: Proceeding with the {len(success)} successful ligands.")
