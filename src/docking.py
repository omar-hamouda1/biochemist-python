"""
docking.py — Utilities for standardized molecular docking with Smina.

This module provides helper functions to parse standardized Smina docking
outputs and rank docking results for the virtual screening pipeline.
"""

import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_smina_output(pdbqt_path: str) -> List[float]:
    """Parse standardized Smina PDBQT output.

    Reads ``REMARK minimizedAffinity`` lines from a docked PDBQT file
    and returns the reported docking affinity scores (kcal/mol).

    Parameters
    ----------
    pdbqt_path : str
        Path to the docked ``*_out.pdbqt`` file.

    Returns
    -------
    list of float
        Docking affinity scores in kcal/mol, with more negative values
        representing better predicted docking scores.
        Returns an empty list if the file is missing, empty, or cannot
        be parsed.
    """
    affinities = []
    path = Path(pdbqt_path)

    if not path.exists() or path.stat().st_size == 0:
        return affinities

    with path.open() as handle:
        for line in handle:
            if line.startswith("REMARK minimizedAffinity"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        affinities.append(float(parts[2]))
                    except ValueError:
                        continue

    return affinities


def load_affinities(
    csv_path: str = "docking/results/standardized_affinities.csv",
) -> pd.DataFrame:
    """Load the standardized docking affinity report.

    Parameters
    ----------
    csv_path : str
        Path to the standardized docking CSV report.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the standardized docking results.
    """
    return pd.read_csv(csv_path)


def get_top_hits(
    csv_path: str = "docking/results/standardized_affinities.csv",
    n: int = 10,
) -> pd.DataFrame:
    """Return the top N ligands ranked by docking affinity.

    Parameters
    ----------
    csv_path : str
        Path to the standardized docking CSV report.
    n : int, default 10
        Number of top hits to return.

    Returns
    -------
    pandas.DataFrame
        Top N rows sorted by affinity, with the most negative score first.
    """
    df = load_affinities(csv_path)
    return df.sort_values("affinity").head(n)
