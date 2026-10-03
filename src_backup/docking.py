"""
docking.py — Utilities for molecular docking with AutoDock Vina / Smina.

This module provides helper functions to prepare, run, and parse
docking results for the virtual screening pipeline.

Example
-------
>>> from src.docking import parse_vina_output, get_top_hits
>>> affinities = parse_vina_output("docking/results/13U_out.pdbqt")
>>> print(f"Best affinity: {affinities[0]:.2f} kcal/mol")
"""

import os
import re
import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Tuple

import pandas as pd


def parse_vina_output(pdbqt_path: str) -> List[float]:
    """Parse AutoDock Vina PDBQT output and extract binding affinities.

    Reads the ``REMARK VINA RESULT`` lines from a docked PDBQT file
    and returns a list of affinities (kcal/mol) for each pose.

    Parameters
    ----------
    pdbqt_path : str
        Path to the docked ``*_out.pdbqt`` file.

    Returns
    -------
    list of float
        Binding affinities in kcal/mol (negative = better binding).
        Empty list if file is empty or cannot be parsed.

    Example
    -------
    >>> affinities = parse_vina_output("docking/results/13U_out.pdbqt")
    >>> print(f"Best pose: {affinities[0]:.1f} kcal/mol")
    Best pose: -9.5 kcal/mol
    """
    affinities = []
    path = Path(pdbqt_path)

    if not path.exists() or path.stat().st_size == 0:
        return affinities

    with open(path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.strip().split()
                if len(parts) >= 4:
                    try:
                        affinities.append(float(parts[3]))
                    except ValueError:
                        continue
    return affinities


def load_affinities(csv_path: str = "docking/results/affinities_full.csv") -> pd.DataFrame:
    """Load the full docking affinities CSV.

    Parameters
    ----------
    csv_path : str
        Path to the affinities CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ``rank``, ``ligand_id``, ``affinity``.
    """
    return pd.read_csv(csv_path)


def get_top_hits(
    csv_path: str = "docking/results/affinities_full.csv",
    n: int = 10,
) -> pd.DataFrame:
    """Return the top N hits ranked by binding affinity.

    Parameters
    ----------
    csv_path : str
        Path to the affinities CSV file.
    n : int, default 10
        Number of top hits to return.

    Returns
    -------
    pd.DataFrame
        Top N rows sorted by affinity (most negative first).
    """
    df = load_affinities(csv_path)
    return df.sort_values("affinity").head(n)


def read_box_config(config_path: str = "docking/box_config.txt") -> Dict[str, float]:
    """Read docking box configuration from a text file.

    Parameters
    ----------
    config_path : str
        Path to the box configuration file.

    Returns
    -------
    dict
        Dictionary with keys: ``center_x``, ``center_y``, ``center_z``,
        ``size_x``, ``size_y``, ``size_z``.

    Example
    -------
    >>> box = read_box_config()
    >>> print(f"Center: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
    """
    config = {}
    with open(config_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            match = re.match(r"(\w+)\s*=\s*([\d.\-]+)", line)
            if match:
                config[match.group(1)] = float(match.group(2))
    return config
