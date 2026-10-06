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
    """Convert the standardized PDBQT pose directly to SDF."""
    result = subprocess.run(
        ["obabel", "-ipdbqt", pdbqt_path, "-osdf", "-O", output_sdf],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0 or not os.path.exists(output_sdf):
        detail = (result.stderr or result.stdout).strip()
        return False, f"Open Babel failed: {detail[:160]}"

    return True, "PDBQT -> SDF"


def make_graph_agnostic_mol(mol):
    """
    Make a comparison copy that ignores formal charge and bond order while
    retaining atom elements and connectivity.
    """
    out = Chem.Mol(mol)

    for atom in out.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(True)

    editable = Chem.RWMol(out)
    for bond in editable.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)

    return editable.GetMol()


def rebuild_from_template(template_path, pose_path, output_sdf):
    """
    Reconstruct the ligand using the original SDF topology/charges and the
    docked PDBQT-derived coordinates.
    """
    template = Chem.MolFromMolFile(
        template_path,
        removeHs=False,
        sanitize=False,
    )
    pose = Chem.MolFromMolFile(
        pose_path,
        removeHs=False,
        sanitize=False,
    )

    if template is None:
        return False, "Cannot read template"
    if pose is None:
        return False, "Cannot read converted pose"

    template_hfree = remove_all_hydrogens(template)
    pose_hfree = remove_all_hydrogens(pose)

    if template_hfree.GetNumAtoms() != pose_hfree.GetNumAtoms():
        return (
            False,
            f"Heavy-atom count mismatch "
            f"({pose_hfree.GetNumAtoms()} vs {template_hfree.GetNumAtoms()})",
        )

    template_cmp = make_graph_agnostic_mol(template_hfree)
    pose_cmp = make_graph_agnostic_mol(pose_hfree)

    mapping = pose_cmp.GetSubstructMatch(template_cmp)

    if not mapping:
        return False, "No graph isomorphism between template and pose"

    template_symbols = [a.GetSymbol() for a in template_hfree.GetAtoms()]
    pose_symbols = [a.GetSymbol() for a in pose_hfree.GetAtoms()]

    for template_idx, pose_idx in enumerate(mapping):
        if template_symbols[template_idx] != pose_symbols[pose_idx]:
            return False, "Graph mapping produced element mismatch"

    final_mol = Chem.Mol(template)
    if final_mol.GetNumConformers():
        final_conf = final_mol.GetConformer()
    else:
        final_conf = Chem.Conformer(final_mol.GetNumAtoms())
        final_mol.AddConformer(final_conf, assignId=True)

    pose_conf = pose.GetConformer()
    pose_hfree_conf = pose_hfree.GetConformer()
    template_h_to_pose_h = {}

    for template_heavy_idx, pose_heavy_idx in enumerate(mapping):
        pos = pose_hfree_conf.GetAtomPosition(pose_heavy_idx)
        final_conf.SetAtomPosition(template_heavy_idx, pos)

        template_atom = template.GetAtomWithIdx(template_heavy_idx)
        pose_atom = pose.GetAtomWithIdx(pose_heavy_idx)

        template_h = [
            n.GetIdx()
            for n in template_atom.GetNeighbors()
            if n.GetSymbol() == "H"
        ]
        pose_h = [
            n.GetIdx()
            for n in pose_atom.GetNeighbors()
            if n.GetSymbol() == "H"
        ]

        for t_idx, p_idx in zip(template_h, pose_h):
            template_h_to_pose_h[t_idx] = p_idx

    for template_h_idx, pose_h_idx in template_h_to_pose_h.items():
        pos = pose_conf.GetAtomPosition(pose_h_idx)
        final_conf.SetAtomPosition(template_h_idx, pos)

    try:
        Chem.SanitizeMol(final_mol)
    except Exception as exc:
        return False, f"Template reconstruction sanitization failed: {exc}"

    writer = Chem.SDWriter(output_sdf)
    writer.write(final_mol)
    writer.close()

    return True, "template topology + docked coordinates"


def validate_pose(ligand_id, pose_sdf):
    """Validate the final SDF with RDKit."""
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
        return False, "Cannot read final pose"

    if element_counts(template) != element_counts(pose):
        return False, "Element composition mismatch"

    try:
        Chem.SanitizeMol(pose)
    except Exception as exc:
        return False, f"RDKit sanitization failed: {exc}"

    return True, "validated"


success, failed = [], []

print("Preparing standardized Top 5 for ProLIF...\n")

for lid in top_5:
    pdbqt = f"{POSE_DIR}/{lid}_out.pdbqt"
    output = f"{OUTPUT_DIR}/{lid}_fixed.sdf"

    if not os.path.exists(pdbqt):
        print(f"  FAIL {lid:5s} — standardized pose not found")
        failed.append(lid)
        continue

    template = f"ligands/{lid}.sdf"

    with tempfile.NamedTemporaryFile(
        suffix=".sdf",
        dir=OUTPUT_DIR,
        delete=False,
    ) as tmp:
        converted = tmp.name

    ok, method = convert_pdbqt_to_sdf(pdbqt, converted)

    if not ok:
        print(f"  FAIL {lid:5s} — {method}")
        failed.append(lid)
        try:
            os.remove(converted)
        except OSError:
            pass
        continue

    ok, validation = validate_pose(lid, converted)

    if ok:
        os.replace(converted, output)
        print(f"  OK   {lid:5s} via {method}; {validation}")
        success.append(lid)
        continue

    ok, rebuild_method = rebuild_from_template(
        template, converted, output
    )

    try:
        os.remove(converted)
    except OSError:
        pass

    if ok:
        print(f"  OK   {lid:5s} via {rebuild_method}; validated")
        success.append(lid)
    else:
        print(f"  FAIL {lid:5s} — {validation}; fallback: {rebuild_method}")
        failed.append(lid)

print(f"\n{'=' * 60}")
print(f"  Success: {len(success)}/{len(top_5)}")
print(f"  Failed:  {len(failed)}/{len(top_5)}")
print(f"{'=' * 60}")

if failed:
    print(f"\nFailed ligands: {failed}")
