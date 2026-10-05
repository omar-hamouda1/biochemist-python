import os

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


def pdbqt_element(atom_name, atom_type):
    """Map common PDBQT atom types/names to chemical elements."""
    atom_type = atom_type.strip().upper()
    atom_name = atom_name.strip().upper()

    if atom_type.startswith("H") or atom_name.startswith("H"):
        return "H"

    type_map = {
        "A": "C",
        "C": "C",
        "N": "N",
        "NA": "N",
        "OA": "O",
        "O": "O",
        "SA": "S",
        "S": "S",
        "P": "P",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
    }
    if atom_type in type_map:
        return type_map[atom_type]

    if atom_name.startswith("CL"):
        return "Cl"
    if atom_name.startswith("BR"):
        return "Br"

    return atom_name[:1].upper()


def parse_pdbqt_heavy_atoms(pdbqt_path):
    """
    Read heavy-atom coordinates directly from PDBQT.

    Hydrogen detection uses PDBQT atom type/name rather than the element
    columns because PDBQT stores AutoDock atom types such as A and NA.
    """
    atoms = []

    with open(pdbqt_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            atom_type = line[77:79].strip()

            element = pdbqt_element(atom_name, atom_type)
            if element == "H":
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atoms.append(
                {
                    "name": atom_name,
                    "type": atom_type,
                    "element": element,
                    "x": x,
                    "y": y,
                    "z": z,
                }
            )

    return atoms


def build_from_template_and_pose(pdbqt_path, template_sdf, output_sdf):
    """
    Keep chemical connectivity/bond orders from the trusted SDF template and
    transfer docked heavy-atom coordinates from the standardized PDBQT pose.
    """
    try:
        template = Chem.MolFromMolFile(
            template_sdf,
            removeHs=False,
            sanitize=False,
        )
        if template is None:
            return False, "Cannot read template"

        template_noH = remove_all_hydrogens(template)
        pose_atoms = parse_pdbqt_heavy_atoms(pdbqt_path)

        expected = template_noH.GetNumAtoms()
        observed = len(pose_atoms)

        if observed != expected:
            return (
                False,
                f"Heavy-atom mismatch ({observed} vs {expected})",
            )

        template_symbols = [
            atom.GetSymbol()
            for atom in template_noH.GetAtoms()
        ]
        pose_symbols = [atom["element"] for atom in pose_atoms]

        if template_symbols != pose_symbols:
            return (
                False,
                "Heavy-atom element/order mismatch "
                f"(template={template_symbols}, pose={pose_symbols})",
            )

        final_mol = Chem.Mol(template_noH)
        final_mol.RemoveAllConformers()

        conf = Chem.Conformer(expected)
        for idx, atom in enumerate(pose_atoms):
            conf.SetAtomPosition(
                idx,
                (atom["x"], atom["y"], atom["z"]),
            )

        final_mol.AddConformer(conf, assignId=True)
        final_mol = Chem.AddHs(final_mol, addCoords=True)
        Chem.SanitizeMol(final_mol)

        writer = Chem.SDWriter(output_sdf)
        writer.write(final_mol)
        writer.close()

        return True, "template + PDBQT coordinates"

    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


success, failed = [], []
print("Preparing standardized Top 5 for ProLIF...\n")

for lid in top_5:
    pdbqt = f"{POSE_DIR}/{lid}_out.pdbqt"
    output = f"{OUTPUT_DIR}/{lid}_fixed.sd"

    if not os.path.exists(pdbqt):
        print(f"  FAIL {lid:5s} — standardized pose not found")
        failed.append(lid)
        continue

    template = f"ligands/{lid}.sdf"

    ok, method = build_from_template_and_pose(
        pdbqt, template, output
    )

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
