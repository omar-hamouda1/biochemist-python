"""Cross-file consistency checks for the current standardized pipeline."""

import ast
import csv
import json
from pathlib import Path

import pandas as pd

from src.docking_config import load_docking_config


ROOT = Path(__file__).resolve().parents[1]


def read_ids(path: Path):
    with path.open(newline="") as handle:
        return [
            row["ligand_id"].strip()
            for row in csv.DictReader(handle)
            if row.get("ligand_id")
        ]


def test_current_result_chain_is_coherent():
    config = load_docking_config()

    candidates = {path.stem for path in (ROOT / "ligands").glob("*.sdf")}
    validated = set(read_ids(config.manifest))
    exceptions = set(read_ids(config.report.parent / "docking_exceptions.csv"))

    assert len(candidates) == 117
    assert len(validated) == config.expected_ligands == 111
    assert len(exceptions) == 6
    assert candidates == validated | exceptions
    assert not validated & exceptions

    docking = pd.read_csv(config.report)
    assert len(docking) == 111
    assert docking["ligand_id"].is_unique
    assert set(docking["status"]).issubset({"ok", "existing"})
    assert docking["status"].notna().all()
    assert set(docking["ligand_id"]) == validated
    assert docking["affinity"].notna().all()

    docking = docking.sort_values(["affinity", "ligand_id"], kind="mergesort").reset_index(drop=True)
    top5 = docking.head(5)["ligand_id"].tolist()
    assert top5 == ["R11", "13U", "BAH", "12U", "607"]

    prolif = pd.read_csv(ROOT / "docking/results/prolif_summary.csv")
    admet = pd.read_csv(ROOT / "docking/results/admet/admet_top5.csv")
    final = pd.read_csv(ROOT / "results/top_hits_summary.csv")

    assert set(prolif["ligand_id"]) == set(top5)
    assert set(admet["ligand_id"]) == set(top5)

    final_top5 = final.sort_values("rank").head(5)
    assert final["ligand_id"].is_unique
    assert len(final) == 111
    assert final_top5["ligand_id"].tolist() == top5

    final_scores = dict(zip(final["ligand_id"], final["affinity"]))
    docking_scores = dict(zip(docking["ligand_id"], docking["affinity"]))
    for ligand_id in top5:
        assert final_scores[ligand_id] == docking_scores[ligand_id]


def test_current_analysis_notebooks_have_valid_python_code():
    for notebook_path in [
        ROOT / "notebooks/14_molecular_docking.ipynb",
        ROOT / "notebooks/solutions/14b_docking_solved.ipynb",
    ]:
        with notebook_path.open() as handle:
            notebook = json.load(handle)

        for index, cell in enumerate(notebook["cells"], start=1):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            ast.parse(source, filename=f"{notebook_path}:{index}")

def test_current_top5_artifact_contract_is_sdf():
    """Downstream Top-5 scripts and tracked artifacts must use one extension."""
    top5 = ["R11", "13U", "BAH", "12U", "607"]
    fixed_dir = ROOT / "docking/results/fixed"

    for ligand_id in top5:
        assert (fixed_dir / f"{ligand_id}_fixed.sdf").exists()
        assert not (fixed_dir / f"{ligand_id}_fixed.sd").exists()

    for script_name in [
        "scripts/prepare_top5.py",
        "scripts/prolif_analysis.py",
        "scripts/admet_top5.py",
    ]:
        source = (ROOT / script_name).read_text()
        assert "_fixed.sdf" in source
        assert "_fixed.sd" not in source

