"""
Biochemist Python — Computational Drug Discovery Pipeline.

src package providing reusable modules for:
- Protein structure handling (protein.py)
- ADMET property prediction (admet.py)
- Molecular docking utilities (docking.py)
- 3D visualization helpers (visualization.py)
"""

from src.protein import load_protein, clean_protein, get_residue_count, extract_chain
from src.admet import calculate_lipinski_rules
from src.docking import parse_vina_output, get_top_hits, load_affinities, read_box_config

__all__ = [
    # Protein
    "load_protein",
    "clean_protein",
    "get_residue_count",
    "extract_chain",
    # ADMET
    "calculate_lipinski_rules",
    # Docking
    "parse_vina_output",
    "get_top_hits",
    "load_affinities",
    "read_box_config",
]

__version__ = "1.0.0"
