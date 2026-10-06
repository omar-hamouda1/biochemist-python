"""Prepare the authoritative Trypsin receptor for standardized docking."""

from __future__ import annotations

import subprocess
from pathlib import Path

from src.docking_config import load_docking_config
from src.provenance import (
    build_artifact_metadata,
    executable_version,
    provenance_path,
    write_metadata,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT_ROOT / "pdb" / "protein_h.pdb"
OUTPUT = PROJECT_ROOT / "docking" / "receptor" / "2zq2_receptor.pdbqt"
SCRIPT_PATH = Path(__file__).resolve()
PDB2PQR_NOTEBOOK = PROJECT_ROOT / "notebooks" / "13_binding_site.ipynb"


def is_valid_pdbqt(path: Path) -> bool:
    """Return True when the receptor PDBQT contains atom records."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open(encoding="utf-8", errors="replace") as handle:
        return any(line.startswith(("ATOM", "HETATM")) for line in handle)


def prepare_receptor() -> None:
    """Convert the prepared PDB2PQR protein to a rigid receptor PDBQT."""
    load_docking_config()

    if not SOURCE.exists():
        raise FileNotFoundError(f"Missing prepared protein: {SOURCE}")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.unlink(missing_ok=True)
    provenance_path(OUTPUT).unlink(missing_ok=True)

    command = [
        "obabel",
        "-ipdb",
        str(SOURCE),
        "-opdbqt",
        "-O",
        str(OUTPUT),
        "-xr",
    ]
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not is_valid_pdbqt(OUTPUT):
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"Open Babel receptor conversion failed: {detail[:500]}")

    metadata = build_artifact_metadata(
        source_path=SOURCE,
        producer_path=SCRIPT_PATH,
        tool_versions={
            "openbabel": executable_version("obabel"),
            "pdb2pqr": executable_version("pdb2pqr"),
        },
        protocol={
            "input": str(SOURCE.relative_to(PROJECT_ROOT)),
            "output": str(OUTPUT.relative_to(PROJECT_ROOT)),
            "command": command,
            "rigid_receptor": True,
            "protein_preparation_notebook": str(
                PDB2PQR_NOTEBOOK.relative_to(PROJECT_ROOT)
            ),
        },
        artifact_kind="prepared_receptor_pdbqt",
    )
    write_metadata(
        provenance_path(OUTPUT),
        metadata,
        output_path=OUTPUT,
        extra={
            "upstream_producer": str(PDB2PQR_NOTEBOOK.relative_to(PROJECT_ROOT)),
        },
    )

    print(f"Prepared receptor: {OUTPUT}")
    print(f"Provenance: {provenance_path(OUTPUT)}")


if __name__ == "__main__":
    prepare_receptor()
