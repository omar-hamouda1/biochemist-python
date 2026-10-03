# 🔬 Pipeline Architecture — Detailed Walkthrough

This document explains the full computational drug discovery pipeline
implemented in this project. Each stage maps to one or more Jupyter notebooks.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    COMPUTATIONAL DRUG DISCOVERY                     │
│                        Target: Trypsin (2ZQ2)                       │
└─────────────────────────────────────────────────────────────────────┘

Stage 1: Foundations          Stage 2: Structure          Stage 3: Discovery
┌──────────────────┐         ┌──────────────────┐       ┌──────────────────┐
│ Python Basics    │         │ Biopython/mmCIF  │       │ EC Class Search  │
│ File Parsing     │────────▶│ RCSB Web API     │──────▶│ Ligand Library   │
│ Data Analysis    │         │ 3D Visualization │       │ (117 ligands)    │
│ Regression       │         │ Structure QC     │       │ Binding Site     │
└──────────────────┘         └──────────────────┘       └──────────────────┘
 NB: 01-08                    NB: 09-11                   NB: 12-13
                                                              │
                                                              ▼
Stage 6: Analysis            Stage 5: ADMET             Stage 4: Docking
┌──────────────────┐         ┌──────────────────┐       ┌──────────────────┐
│ MD Trajectory    │         │ Lipinski Ro5     │       │ AutoDock Vina    │
│ RMSD / RMSF      │◀────────│ Druglikeness     │◀──────│ Virtual Screen   │
│ MMPBSA           │         │ ADMET Filters    │       │ 117 ligands      │
│ ProLIF Heatmap   │         │                  │       │ Ranked by ΔG     │
└──────────────────┘         └──────────────────┘       └──────────────────┘
 NB: 16                       NB: 15                      NB: 14
```

---

## Stage 1: Python Foundations (Notebooks 01–08)

### What You Learn
- Python syntax, variables, loops, functions, conditionals
- File I/O: reading CSV and PDB files
- Processing multiple files with `glob` and `os`
- Data manipulation with `pandas`
- Linear regression (Bradford protein assay)
- Plotting with `matplotlib` and `seaborn`
- Nonlinear regression (Michaelis-Menten enzyme kinetics)

### Key Concepts
| Notebook | Biochemistry Concept | Python Concept |
|----------|---------------------|----------------|
| 01 | — | Variables, types, loops |
| 02 | PDB file format | String parsing, file I/O |
| 03 | Batch analysis | `glob`, `os.path` |
| 04 | Data tables | `pandas` DataFrames |
| 05 | Bradford assay | `scipy.stats.linregress` |
| 06 | Assay visualization | `matplotlib`, `seaborn` |
| 07 | Michaelis-Menten | `scipy.optimize.curve_fit` |
| 08 | Inhibition kinetics | Nonlinear regression |

---

## Stage 2: Structural Bioinformatics (Notebooks 09–11)

### What You Learn
- Parsing mmCIF files with Biopython
- Querying the RCSB PDB via REST API
- 3D structure visualization with ICN3D and py3Dmol
- Understanding resolution, R-factor, and structure quality

### Key Data
- **Myoglobin structures**: 40+ CIF files in `pdb_files/`
- **Target protein**: Trypsin (2ZQ2) — 1.7 Å resolution

---

## Stage 3: Drug Discovery (Notebooks 12–13)

### What You Learn
- Searching PDB by EC classification (EC 3.4.21 → Serine Proteases)
- Extracting bound ligands from crystal structures
- Building a screening library (117 unique ligands)
- Binding site identification and visualization
- Understanding protein-ligand interactions

### Output
- `ligands/` directory: 117 SDF files
- Binding site residue list
- Reference ligand (13U) interaction profile

---

## Stage 4: Virtual Screening (Notebook 14)

### What You Learn
- Receptor preparation (PDB → PDBQT)
- Ligand preparation (SDF → PDBQT)
- Docking box definition (centered on binding site)
- Running AutoDock Vina / Smina
- Parsing and ranking results

### Configuration
```
Box center: (17.672, -8.256, 10.688) Å
Box size: 25.0 × 25.0 × 25.0 Å
Exhaustiveness: 8
```

### Top 5 Results
| Rank | Ligand | Affinity (kcal/mol) |
|------|--------|---------------------|
| 1 | 13U | -9.52 |
| 2 | R11 | -9.51 |
| 3 | BAH | -9.40 |
| 4 | 607 | -9.40 |
| 5 | 12U | -9.35 |

---

## Stage 5: ADMET Prediction (Notebook 15)

### What You Learn
- Lipinski's Rule of Five
- Molecular descriptors (MW, LogP, HBD, HBA)
- Drug-likeness assessment
- Filtering hits by ADMET properties

---

## Stage 6: Molecular Dynamics (Notebook 16)

### What You Learn
- System preparation (protein + ligand + solvent + ions)
- Energy minimization
- NVT and NPT equilibration
- Production MD run
- Trajectory analysis (RMSD, RMSF, hydrogen bonds)
- MMPBSA binding free energy decomposition
- ProLIF interaction fingerprint heatmaps

### Software
- **OpenMM** for simulation
- **MDAnalysis** for trajectory analysis
- **ProLIF** for interaction fingerprints

---

## Directory Map

```
biochemist-python_ORGANIZED/
├── notebooks/          # 16 Jupyter notebooks (main pipeline)
├── src/                # Reusable Python modules
│   ├── protein.py      #   Protein structure utilities
│   ├── admet.py        #   ADMET/Lipinski calculations
│   ├── docking.py      #   Docking result parsing
│   └── visualization.py#   3D visualization helpers
├── scripts/            # Standalone analysis scripts
├── tests/              # Unit tests (pytest)
├── data/               # Raw experimental data (CSV)
├── pdb/                # Protein structure files
├── pdb_files/          # Myoglobin CIF collection
├── ligands/            # 117 ligand SDF files
├── docking/            # Docking inputs and outputs
├── md/                 # MD simulation files
├── figures/            # Generated publication figures
├── results/            # Final merged results
├── configs/            # Configuration files
├── molssi_data/        # MolSSI workshop reference
└── docs/               # Documentation
```
