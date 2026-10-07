"""Prepare the Trypsin protein for PDB2PQR and downstream docking.

This script is the reusable, authoritative protein-preparation entry point.
It reproduces the structure-processing logic introduced in Notebook 13 while
avoiding non-standard PDB header/REMARK records that can trigger PDB2PQR
parsing errors.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import MDAnalysis as mda


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT_ROOT / "pdb" / "2zq2.pdb"
RESOLVED = PROJECT_ROOT / "pdb" / "protein_2zq2.pdb"
PQR = PROJECT_ROOT / "pdb" / "protein.pqr"
PREPARED = PROJECT_ROOT / "pdb" / "protein_h.pdb"


def select_dominant_protein_conformer(protein_group):
    """Select the highest-mean-occupancy conformer for each protein residue."""
    selected_indices = []
    decisions = []

    for residue in protein_group.residues:
        altlocs = sorted(
            {
                atom.altLoc.strip()
                for atom in residue.atoms
                if atom.altLoc and atom.altLoc.strip()
            }
        )

        if not altlocs:
            selected_indices.extend(int(atom.index) for atom in residue.atoms)
            continue

        mean_occupancy = {}
        for altloc in altlocs:
            occupancies = [
                float(atom.occupancy)
                for atom in residue.atoms
                if atom.altLoc.strip() == altloc
                and atom.occupancy is not None
            ]
            if not occupancies:
                raise ValueError(
                    f"No occupancy values found for "
                    f"{residue.resname} {residue.resid} altLoc {altloc}"
                )
            mean_occupancy[altloc] = sum(occupancies) / len(occupancies)

        selected_altloc = max(
            mean_occupancy,
            key=lambda alt: (mean_occupancy[alt], alt == "A"),
        )

        for atom in residue.atoms:
            atom_altloc = atom.altLoc.strip()
            if not atom_altloc or atom_altloc == selected_altloc:
                selected_indices.append(int(atom.index))

        decisions.append(
            {
                "resname": residue.resname,
                "resid": str(residue.resid),
                "selected": selected_altloc,
                "occupancy": mean_occupancy,
            }
        )

    selected_indices = sorted(set(selected_indices))
    return protein_group.universe.atoms[selected_indices], decisions


def write_pdb2pqr_safe_pdb(atom_group, output_path: Path) -> None:
    """Write coordinate records only, with altLoc labels cleared."""
    temporary_path = output_path.with_suffix(".tmp.pdb")
    atom_group.write(str(temporary_path))

    written = 0
    last_line = ""
    with temporary_path.open(encoding="utf-8") as source, output_path.open(
        "w", encoding="utf-8"
    ) as target:
        for line in source:
            last_line = line
            if line.startswith(("ATOM  ", "HETATM")):
                if len(line) >= 17 and line[16] != " ":
                    line = line[:16] + " " + line[17:]
                target.write(line)
                written += 1
            elif line.startswith(("TER   ", "END   ")):
                target.write(line)

        if not last_line.startswith("END"):
            target.write("END\n")

    temporary_path.unlink(missing_ok=True)

    if written != atom_group.n_atoms:
        raise RuntimeError(
            f"Coordinate-record count mismatch: wrote {written}, "
            f"expected {atom_group.n_atoms}"
        )


def run_pdb2pqr() -> None:
    """Run PDB2PQR with the project's fixed preparation protocol."""
    command = [
        "pdb2pqr",
        "--ff",
        "PARSE",
        f"--pdb-output={PREPARED}",
        "--pH=7.4",
        str(RESOLVED),
        str(PQR),
    ]
    result = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    combined = (result.stdout or "") + ("\n" + result.stderr if result.stderr else "")
    if result.returncode != 0:
        raise RuntimeError(
            "PDB2PQR failed:\n" + combined[-4000:]
        )

    if not PQR.exists() or PQR.stat().st_size == 0:
        raise RuntimeError("PDB2PQR produced an empty or missing PQR file")
    if not PREPARED.exists() or PREPARED.stat().st_size == 0:
        raise RuntimeError("PDB2PQR produced an empty or missing prepared PDB")


def prepare_protein() -> None:
    """Resolve altLocs, write a clean PDB, and run PDB2PQR."""
    if not SOURCE.exists():
        raise FileNotFoundError(f"Missing source PDB: {SOURCE}")

    universe = mda.Universe(str(SOURCE))
    protein = universe.select_atoms("protein")
    if protein.n_atoms == 0:
        raise RuntimeError("Source PDB contains no protein atoms")

    protein_dominant, decisions = select_dominant_protein_conformer(protein)
    RESOLVED.parent.mkdir(parents=True, exist_ok=True)

    write_pdb2pqr_safe_pdb(protein_dominant, RESOLVED)

    print("=" * 70)
    print("  Protein preparation: Trypsin 2ZQ2")
    print("=" * 70)
    print(f"Source atoms:   {protein.n_atoms:,}")
    print(f"Resolved atoms: {protein_dominant.n_atoms:,}")
    for decision in decisions:
        occupancy_text = ", ".join(
            f"{alt}={value:.2f}"
            for alt, value in sorted(decision["occupancy"].items())
        )
        print(
            f"  {decision['resname']} {decision['resid']}: "
            f"{occupancy_text} -> selected {decision['selected']}"
        )

    print(f"\nWrote clean PDB: {RESOLVED}")
    print("Coordinate records only; non-standard REMARK/header records are excluded.")

    run_pdb2pqr()

    print(f"Wrote PQR:        {PQR}")
    print(f"Wrote prepared PDB:{PREPARED}")


if __name__ == "__main__":
    prepare_protein()
