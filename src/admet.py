"""
admet.py — ADMET property prediction using RDKit descriptors.

Calculates Lipinski's Rule of Five and additional druglikeness properties
for candidate molecules identified through virtual screening.

Example
-------
>>> from src.admet import calculate_lipinski_rules
>>> df = calculate_lipinski_rules(["CC(=O)Oc1ccccc1C(=O)O"])  # Aspirin
>>> print(df[["MW", "LogP", "Passes_Ro5"]])
"""

from typing import List, Dict, Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def calculate_lipinski_rules(smiles_list: List[str]) -> pd.DataFrame:
    """Calculate Lipinski's Rule of Five properties for a list of SMILES.

    Lipinski's Rule of Five predicts oral bioavailability:
    - Molecular Weight ≤ 500 Da
    - LogP ≤ 5
    - H-bond Donors ≤ 5
    - H-bond Acceptors ≤ 10

    Parameters
    ----------
    smiles_list : list of str
        SMILES strings of molecules to evaluate.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: ``SMILES``, ``Valid``, ``MW``, ``LogP``,
        ``HBD``, ``HBA``, ``Passes_Ro5``.
        Invalid SMILES will have ``Valid=False`` and NaN for properties.

    Example
    -------
    >>> df = calculate_lipinski_rules(["CC(=O)Oc1ccccc1C(=O)O"])
    >>> assert df.iloc[0]["Passes_Ro5"] == True
    >>> assert 170 < df.iloc[0]["MW"] < 190  # Aspirin ≈ 180 Da
    """
    results: List[Dict[str, Any]] = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({"SMILES": smiles, "Valid": False})
            continue

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        passes_lipinski = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)

        results.append({
            "SMILES": smiles,
            "Valid": True,
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "Passes_Ro5": passes_lipinski,
        })

    return pd.DataFrame(results)


if __name__ == "__main__":
    # Quick smoke test
    test_smiles = ["CC(=O)Oc1ccccc1C(=O)O"]  # Aspirin
    df = calculate_lipinski_rules(test_smiles)
    print(df)
