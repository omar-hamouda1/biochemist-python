"""
visualization.py — 3D molecular visualization helpers.

Provides convenience wrappers for py3Dmol to visualize proteins,
ligands, and docking poses inside Jupyter notebooks.

Example
-------
>>> from src.visualization import view_protein_ligand
>>> view = view_protein_ligand("pdb/protein_2zq2.pdb", "pdb/ligand_A.pdb")
>>> view.show()
"""

from pathlib import Path
from typing import Optional

try:
    import py3Dmol
except ImportError:
    py3Dmol = None


def _check_py3dmol():
    """Raise ImportError if py3Dmol is not installed."""
    if py3Dmol is None:
        raise ImportError(
            "py3Dmol is required for visualization. "
            "Install it with: pip install py3Dmol"
        )


def view_protein(
    pdb_path: str,
    style: str = "cartoon",
    color: str = "spectrum",
    width: int = 800,
    height: int = 500,
):
    """Render a protein structure in 3D.

    Parameters
    ----------
    pdb_path : str
        Path to the protein PDB file.
    style : str, default "cartoon"
        Representation style (``cartoon``, ``stick``, ``surface``, ``sphere``).
    color : str, default "spectrum"
        Color scheme (``spectrum``, ``chain``, ``ssType``).
    width : int, default 800
        Viewer width in pixels.
    height : int, default 500
        Viewer height in pixels.

    Returns
    -------
    py3Dmol.view
        The 3D viewer object. Call ``.show()`` to display.
    """
    _check_py3dmol()

    with open(pdb_path, "r") as f:
        pdb_data = f.read()

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")
    view.setStyle({style: {"color": color}})
    view.zoomTo()
    return view


def view_protein_ligand(
    protein_path: str,
    ligand_path: str,
    protein_style: str = "cartoon",
    ligand_style: str = "stick",
    width: int = 800,
    height: int = 500,
):
    """Render a protein-ligand complex in 3D.

    The protein is shown as a cartoon (colored by spectrum) and
    the ligand as green sticks, zoomed to the binding site.

    Parameters
    ----------
    protein_path : str
        Path to the protein PDB file.
    ligand_path : str
        Path to the ligand PDB/SDF file.
    protein_style : str, default "cartoon"
        Representation for the protein.
    ligand_style : str, default "stick"
        Representation for the ligand.
    width : int, default 800
        Viewer width in pixels.
    height : int, default 500
        Viewer height in pixels.

    Returns
    -------
    py3Dmol.view
        The 3D viewer object. Call ``.show()`` to display.

    Example
    -------
    >>> view = view_protein_ligand("pdb/protein_2zq2.pdb", "pdb/ligand_A.pdb")
    >>> view.show()
    """
    _check_py3dmol()

    with open(protein_path, "r") as f:
        prot_data = f.read()
    with open(ligand_path, "r") as f:
        lig_data = f.read()

    # Detect ligand format from extension
    lig_format = "sdf" if ligand_path.endswith(".sdf") else "pdb"

    view = py3Dmol.view(width=width, height=height)
    view.addModel(prot_data, "pdb")
    view.setStyle({"model": 0}, {protein_style: {"color": "spectrum"}})
    view.addModel(lig_data, lig_format)
    view.setStyle(
        {"model": 1},
        {ligand_style: {"colorscheme": "greenCarbon", "radius": 0.2}},
    )
    view.zoomTo({"model": 1})
    return view


def view_binding_site(
    protein_path: str,
    ligand_path: str,
    surface: bool = True,
    width: int = 800,
    height: int = 500,
):
    """Render a binding site with hydrophobic surface and ligand sticks.

    Parameters
    ----------
    protein_path : str
        Path to the protein PDB file.
    ligand_path : str
        Path to the ligand PDB file.
    surface : bool, default True
        If ``True``, show protein surface colored by hydrophobicity.
    width : int, default 800
        Viewer width in pixels.
    height : int, default 500
        Viewer height in pixels.

    Returns
    -------
    py3Dmol.view
        The 3D viewer object.
    """
    _check_py3dmol()

    with open(protein_path, "r") as f:
        prot_data = f.read()
    with open(ligand_path, "r") as f:
        lig_data = f.read()

    view = py3Dmol.view(width=width, height=height)
    view.addModel(prot_data, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum", "opacity": 0.7}})

    if surface:
        view.addSurface(
            py3Dmol.VDW,
            {"opacity": 0.5, "color": "white"},
            {"model": 0},
        )

    lig_format = "sdf" if ligand_path.endswith(".sdf") else "pdb"
    view.addModel(lig_data, lig_format)
    view.setStyle(
        {"model": 1},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.15}},
    )
    view.zoomTo({"model": 1})
    return view
