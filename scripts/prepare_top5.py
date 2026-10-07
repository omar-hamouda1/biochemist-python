"""Prepare the current standardized docking Top 5 for downstream analysis."""

from pathlib import Path
import subprocess
import tempfile

import pandas as pd
from rdkit import Chem

from src.docking_config import load_docking_config

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_ROOT / "docking" / "results" / "fixed"
TOP_N = 5


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


def heavy_atom_indices(mol):
    """Return original atom indices for non-hydrogen atoms."""
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"]


def load_top_hits(report_path: Path, expected_count: int = TOP_N):
    ranking = pd.read_csv(report_path)
    required = {"ligand_id", "affinity", "status"}
    missing = required - set(ranking.columns)
    if missing:
        raise RuntimeError(
            f"Standardized report missing columns: {sorted(missing)}"
        )

    ranking["affinity"] = pd.to_numeric(
        ranking["affinity"], errors="coerce"
    )
    ranking = ranking[
        ranking["status"].isin(["ok", "existing"])
    ].dropna(subset=["affinity"])

    ranking = ranking.sort_values(["affinity", "ligand_id"], kind="mergesort").reset_index(drop=True)

    top = ranking.head(expected_count)
    if len(top) != expected_count:
        raise RuntimeError(
            f"Expected {expected_count} complete Top hits, found {len(top)}"
        )
    if top["ligand_id"].duplicated().any():
        raise RuntimeError("Standardized Top hits contain duplicate ligand IDs")

    return top


def convert_pdbqt_to_sdf(pdbqt_path, output_sdf):
    """Convert the standardized PDBQT pose directly to SDF."""
    result = subprocess.run(
        ["obabel", "-ipdbqt", str(pdbqt_path), "-osdf", "-O", str(output_sdf)],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0 or not output_sdf.exists():
        detail = (result.stderr or result.stdout).strip()
        return False, f"Open Babel failed: {detail[:160]}"

    return True, "PDBQT -> SDF"


def make_graph_agnostic_mol(mol):
    """Create a comparison copy that ignores charge and bond order."""
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
    """Rebuild ligand chemistry from the template with docked coordinates."""
    template = Chem.MolFromMolFile(
        str(template_path), removeHs=False, sanitize=False
    )
    pose = Chem.MolFromMolFile(
        str(pose_path), removeHs=False, sanitize=False
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
            "Heavy-atom count mismatch "
            f"({pose_hfree.GetNumAtoms()} vs {template_hfree.GetNumAtoms()})",
        )

    template_cmp = make_graph_agnostic_mol(template_hfree)
    pose_cmp = make_graph_agnostic_mol(pose_hfree)
    mapping = pose_cmp.GetSubstructMatch(template_cmp)

    if not mapping:
        return False, "No graph isomorphism between template and pose"

    template_symbols = [a.GetSymbol() for a in template_hfree.GetAtoms()]
    pose_symbols = [a.GetSymbol() for a in pose_hfree.GetAtoms()]

    template_heavy_indices = heavy_atom_indices(template)
    pose_heavy_indices = heavy_atom_indices(pose)

    for template_hfree_idx, pose_hfree_idx in enumerate(mapping):
        if template_symbols[template_hfree_idx] != pose_symbols[pose_hfree_idx]:
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

    for template_hfree_idx, pose_hfree_idx in enumerate(mapping):
        template_idx = template_heavy_indices[template_hfree_idx]
        pose_idx = pose_heavy_indices[pose_hfree_idx]

        final_conf.SetAtomPosition(
            template_idx,
            pose_hfree_conf.GetAtomPosition(pose_hfree_idx),
        )

        template_h = [
            n.GetIdx()
            for n in template.GetAtomWithIdx(template_idx).GetNeighbors()
            if n.GetSymbol() == "H"
        ]
        pose_h = [
            n.GetIdx()
            for n in pose.GetAtomWithIdx(pose_idx).GetNeighbors()
            if n.GetSymbol() == "H"
        ]

        for t_idx, p_idx in zip(template_h, pose_h):
            template_h_to_pose_h[t_idx] = p_idx

    for template_h_idx, pose_h_idx in template_h_to_pose_h.items():
        final_conf.SetAtomPosition(
            template_h_idx,
            pose_conf.GetAtomPosition(pose_h_idx),
        )

    try:
        Chem.SanitizeMol(final_mol)
    except Exception as exc:
        return False, f"Template reconstruction sanitization failed: {exc}"

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(final_mol)
    writer.close()

    return True, "template topology + docked coordinates"


def validate_pose(ligand_id, pose_sdf):
    """Validate the final SDF against the source ligand composition."""
    template_path = PROJECT_ROOT / "ligands" / f"{ligand_id}.sdf"

    template = Chem.MolFromMolFile(
        str(template_path), removeHs=False, sanitize=False
    )
    pose = Chem.MolFromMolFile(
        str(pose_sdf), removeHs=False, sanitize=False
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


def prepare_top5():
    """Prepare and validate exactly the current standardized Top 5."""
    config = load_docking_config()
    ranking = load_top_hits(config.report)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    top_5 = ranking["ligand_id"].astype(str).tolist()

    print(f"Standardized Top 5: {top_5}\n")
    success, failed = [], []

    for lid in top_5:
        pdbqt = config.output_dir / f"{lid}_out.pdbqt"
        output = OUTPUT_DIR / f"{lid}_fixed.sdf"
        template = PROJECT_ROOT / "ligands" / f"{lid}.sdf"

        if not pdbqt.exists():
            print(f"  FAIL {lid:5s} — standardized pose not found")
            failed.append(lid)
            continue

        with tempfile.NamedTemporaryFile(
            suffix=".sdf",
            dir=OUTPUT_DIR,
            delete=False,
        ) as tmp:
            converted = Path(tmp.name)

        ok, method = convert_pdbqt_to_sdf(pdbqt, converted)

        if not ok:
            print(f"  FAIL {lid:5s} — {method}")
            failed.append(lid)
            converted.unlink(missing_ok=True)
            continue

        ok, validation = validate_pose(lid, converted)

        if ok:
            converted.replace(output)
            print(f"  OK   {lid:5s} via {method}; {validation}")
            success.append(lid)
            continue

        ok, rebuild_method = rebuild_from_template(
            template, converted, output
        )
        converted.unlink(missing_ok=True)

        if ok:
            print(f"  OK   {lid:5s} via {rebuild_method}; validated")
            success.append(lid)
        else:
            print(
                f"  FAIL {lid:5s} — {validation}; "
                f"fallback: {rebuild_method}"
            )
            failed.append(lid)

    print(f"\n{'=' * 60}")
    print(f"  Success: {len(success)}/{len(top_5)}")
    print(f"  Failed:  {len(failed)}/{len(top_5)}")
    print(f"{'=' * 60}")

    if failed:
        raise RuntimeError(f"Top 5 preparation incomplete: {failed}")


if __name__ == "__main__":
    prepare_top5()
