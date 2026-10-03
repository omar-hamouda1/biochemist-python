"""
protein.py — Utilities for loading, cleaning, and analyzing protein structures.

This module provides helper functions used across the drug discovery pipeline
for handling PDB files with Biopython.

Example
-------
>>> from src.protein import load_protein, clean_protein
>>> structure = load_protein(Path("pdb/2zq2.pdb"))
>>> cleaned = clean_protein(structure, keep_het=False)
"""

from pathlib import Path
from typing import Optional

from Bio import PDB


def load_protein(pdb_path: Path) -> PDB.Structure.Structure:
    """Load a PDB file and return a Biopython Structure object.

    Parameters
    ----------
    pdb_path : Path
        Path to the PDB file (e.g., ``pdb/2zq2.pdb``).

    Returns
    -------
    Bio.PDB.Structure.Structure
        The parsed protein structure.

    Raises
    ------
    FileNotFoundError
        If ``pdb_path`` does not exist.

    Example
    -------
    >>> structure = load_protein(Path("pdb/2zq2.pdb"))
    >>> print(f"Loaded: {structure.id}")
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    return structure


def clean_protein(
    structure: PDB.Structure.Structure,
    keep_het: bool = False,
) -> PDB.Structure.Structure:
    """Remove water, alternate locations, and optionally hetero-atoms.

    This function modifies the structure **in place** and returns it.
    Waters (``HOH``) are always removed.  Hetero-atoms (ligands, ions,
    cofactors) are removed unless ``keep_het=True``.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The protein structure to clean.
    keep_het : bool, default False
        If ``True``, keep hetero-atoms (co-factors, ions, ligands).

    Returns
    -------
    Bio.PDB.Structure.Structure
        The cleaned structure (same object, modified in place).
    """
    for model in structure:
        for chain in model:
            residues_to_remove = []
            for residue in chain:
                het_flag = residue.id[0]

                # Always remove water
                if het_flag == "W":
                    residues_to_remove.append(residue.id)
                # Optionally remove hetero-atoms
                elif not keep_het and het_flag.startswith("H_"):
                    residues_to_remove.append(residue.id)
                else:
                    # Normalize alternate locations — keep only primary conformer
                    for atom in residue:
                        if atom.altloc not in (" ", "A"):
                            atom.altloc = " "

            for res_id in residues_to_remove:
                chain.detach_child(res_id)

    return structure


def get_residue_count(structure: PDB.Structure.Structure) -> int:
    """Count the number of standard amino acid residues in the structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        A Biopython protein structure.

    Returns
    -------
    int
        Number of standard (non-hetero) residues.
    """
    count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    count += 1
    return count


def extract_chain(
    structure: PDB.Structure.Structure,
    chain_id: str = "A",
) -> Optional[PDB.Chain.Chain]:
    """Extract a specific chain from the first model.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The protein structure.
    chain_id : str, default "A"
        The chain identifier to extract.

    Returns
    -------
    Bio.PDB.Chain.Chain or None
        The requested chain, or ``None`` if not found.
    """
    model = structure[0]
    if chain_id in model:
        return model[chain_id]
    return None
