"""
Biochemist Python — Computational Drug Discovery Pipeline.

src package providing reusable modules for:
- Protein structure handling (protein.py)
- ADMET property descriptors (admet.py)
- Standardized molecular docking utilities (docking.py, docking_config.py)
- 3D visualization helpers (visualization.py)
"""

from src.protein import load_protein, clean_protein, get_residue_count, extract_chain
from src.admet import calculate_lipinski_rules
from src.docking import (
    parse_smina_output,
    get_top_hits,
    load_affinities,
)

__all__ = [
    # Protein
    "load_protein",
    "clean_protein",
    "get_residue_count",
    "extract_chain",
    # ADMET
    "calculate_lipinski_rules",
    # Docking
    "parse_smina_output",
    "get_top_hits",
    "load_affinities",
]

__version__ = "1.0.0"
