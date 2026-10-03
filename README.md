# 🧬 Biochemist Python — Computational Drug Discovery Pipeline

> **A complete hands-on workshop:** from Python basics to Molecular Dynamics Simulation.  
> **Target:** Trypsin (`2ZQ2`) — Serine Protease | **Task:** Virtual screening of 117 ligands

[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue?logo=python&logoColor=white)](https://www.python.org/downloads/release/python-3110/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter)](https://jupyter.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2024-green)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey.svg)](LICENSE)
[![CI](https://img.shields.io/badge/CI-GitHub_Actions-blue?logo=githubactions)](https://github.com/features/actions)

---

## 📌 Overview

This project implements a **complete computational drug discovery pipeline** designed for biochemists learning Python. Starting from zero Python knowledge, it walks through every stage of a real drug discovery workflow — all applied to a real protein target.

| | Details |
|---|---|
| 🎯 **Protein Target** | Trypsin — Serine Protease (EC 3.4.21.4) — PDB: [`2ZQ2`](https://www.rcsb.org/structure/2ZQ2) |
| 💊 **Reference Ligand** | `13U` (co-crystallized inhibitor) |
| 📚 **Screening Library** | 117 ligands from PDB (EC class 3.4.21) |
| 🏆 **Top Hit** | `13U` at **-9.52 kcal/mol** (confirmed by MD simulation) |

---

## 🗺️ Pipeline

```
 ┌─────────────────────────────────────────────────────────────────────────┐
 │                                                                         │
 │   Python ──► Data ──► Structures ──► 3D Viz ──► Ligand Discovery       │
 │   Basics     Analysis  (Biopython)   (py3Dmol)   (117 ligands)          │
 │   NB:01-04   NB:05-08  NB:09         NB:10-11    NB:12                  │
 │                                                       │                 │
 │                                                       ▼                 │
 │   Final     MD Sim   ◄── ADMET ◄── Docking ◄── Binding Site           │
 │   Hits      Analysis     Filter     (Vina)      Investigation          │
 │             NB:16        NB:15      NB:14       NB:13                   │
 │                                                                         │
 └─────────────────────────────────────────────────────────────────────────┘
```

> 📖 For a detailed walkthrough of each stage, see [`docs/PIPELINE.md`](docs/PIPELINE.md).

---

## ⚡ Quick Start

```bash
# 1. Clone
git clone https://github.com/<your-username>/biochemist-python.git
cd biochemist-python

# 2. Create environment
conda env create -f environment.yml
conda activate biochem

# 3. Run notebooks
jupyter lab
# Navigate to notebooks/ → start with 01_python_basics.ipynb
```

> 📖 For detailed setup (including WSL, VS Code, and external tools), see [`docs/INSTALLATION.md`](docs/INSTALLATION.md).

---

## 📓 Notebooks

| # | Notebook | Topic | Stage |
|---|----------|-------|-------|
| 01 | [`01_python_basics.ipynb`](notebooks/01_python_basics.ipynb) | Python fundamentals for biochemists | Foundations |
| 02 | [`02_file_parsing.ipynb`](notebooks/02_file_parsing.ipynb) | Parsing CSV and PDB files | Foundations |
| 03 | [`03_multiple_files.ipynb`](notebooks/03_multiple_files.ipynb) | Processing multiple files with `glob` | Foundations |
| 04 | [`04_pandas.ipynb`](notebooks/04_pandas.ipynb) | Data analysis with Pandas | Foundations |
| 05 | [`05_linear_regression.ipynb`](notebooks/05_linear_regression.ipynb) | Bradford protein assay — linear fit | Data Analysis |
| 06 | [`06_plots.ipynb`](notebooks/06_plots.ipynb) | Publication-quality plots | Data Analysis |
| 07 | [`07_nonlinear_regression_part1.ipynb`](notebooks/07_nonlinear_regression_part1.ipynb) | Michaelis-Menten kinetics (Part 1) | Data Analysis |
| 08 | [`08_nonlinear_regression_part2.ipynb`](notebooks/08_nonlinear_regression_part2.ipynb) | Enzyme inhibition kinetics (Part 2) | Data Analysis |
| 09 | [`09_biopython_mmcif.ipynb`](notebooks/09_biopython_mmcif.ipynb) | Structural bioinformatics with Biopython | Structure |
| 10 | [`10_rcsb_web_api.ipynb`](notebooks/10_rcsb_web_api.ipynb) | RCSB PDB REST API queries | Structure |
| 11 | [`11_icn3d_visualization.ipynb`](notebooks/11_icn3d_visualization.ipynb) | 3D structure visualization (iCN3D/py3Dmol) | Structure |
| 12 | [`12_ec_class_ligands.ipynb`](notebooks/12_ec_class_ligands.ipynb) | Ligand discovery by EC class — **117 ligands** | Discovery |
| 13 | [`13_binding_site.ipynb`](notebooks/13_binding_site.ipynb) | Binding site analysis with MDAnalysis | Discovery |
| 14 | [`14_molecular_docking.ipynb`](notebooks/14_molecular_docking.ipynb) | Virtual screening with AutoDock Vina | Docking |
| 15 | [`15_admet_prediction.ipynb`](notebooks/15_admet_prediction.ipynb) | ADMET / Lipinski Ro5 filtering | ADMET |
| 16 | [`16_md_analysis.ipynb`](notebooks/16_md_analysis.ipynb) | MD simulation & trajectory analysis | MD/Analysis |

> **Run order:** Notebooks 01 → 16 sequentially. Each notebook builds on the previous.

---

## 🏆 Key Results

### Top 5 Virtual Screening Hits

| Rank | Ligand | Affinity (kcal/mol) | Passes Ro5 |
|------|--------|---------------------|------------|
| 1 | **13U** | **-9.52** | ✅ |
| 2 | R11 | -9.51 | ✅ |
| 3 | BAH | -9.40 | ✅ |
| 4 | 607 | -9.40 | ✅ |
| 5 | 12U | -9.35 | ✅ |

### MD Simulation Highlights
- **System**: 13U–Trypsin complex in explicit water
- **Key binding residues**: Identified via ProLIF fingerprinting
- **Stability**: Confirmed by RMSD convergence and MMPBSA analysis

---

## 💾 Data Availability

The large MD trajectory files (`*.dcd`, `*.nc`) exceed GitHub's file size limits and are **not** included in this repository.

To run Notebook 16, download the MD simulation data and extract into `md/colab_workshop/02_analysis/`.

---

## 🛠️ Tech Stack

| Category | Tools |
|----------|-------|
| **Core** | Python 3.11, NumPy, Pandas, SciPy |
| **Structure Analysis** | Biopython, MDAnalysis |
| **Cheminformatics** | RDKit, rcsbsearchapi |
| **Docking** | AutoDock Vina / Smina, Open Babel |
| **Interaction Analysis** | ProLIF |
| **Visualization** | py3Dmol, NGLView, iCN3D, Matplotlib, Seaborn |
| **MD Simulation** | OpenMM |
| **ADMET** | RDKit descriptors |
| **Testing** | pytest |

---

## 📁 Project Structure

```
biochemist-python/
├── 📓 notebooks/           # 16 Jupyter notebooks (main pipeline)
│   └── solutions/          # Solved versions of workshop notebooks
├── 📦 src/                 # Reusable Python modules
│   ├── protein.py          #   Protein structure utilities
│   ├── admet.py            #   ADMET / Lipinski calculations
│   ├── docking.py          #   Docking result parsing
│   └── visualization.py    #   3D visualization helpers
├── 🧪 tests/               # Unit tests (pytest)
├── 📜 scripts/             # Standalone analysis scripts
├── 📊 data/                # Raw experimental data (CSV, kinetics)
├── 🧬 pdb/                 # Protein structure files (2ZQ2)
├── 🧬 pdb_files/           # Myoglobin CIF structures (batch)
├── 💊 ligands/             # 117 ligand SDF files
├── 🎯 docking/             # Docking inputs, outputs, configs
├── ⚗️ md/                  # MD simulation (OpenMM)
├── 📈 figures/             # Generated figures
├── 📋 results/             # Final merged results
├── ⚙️ configs/             # Configuration files
├── 📖 docs/                # Documentation
│   ├── PIPELINE.md         #   Detailed pipeline walkthrough
│   └── INSTALLATION.md     #   Setup instructions
├── 📚 molssi_data/         # MolSSI workshop reference
├── environment.yml         # Conda environment specification
├── LICENSE                 # MIT License
├── CONTRIBUTING.md         # Contribution guidelines
├── CODE_OF_CONDUCT.md      # Community standards
└── CHANGELOG.md            # Version history
```

---

## 🧪 Running Tests

```bash
conda activate biochem
pytest tests/ -v
```

---

## 🤝 Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## 📚 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{hamouda2024biochemist,
  author  = {Hamouda, Omar},
  title   = {Biochemist Python: Computational Drug Discovery Pipeline},
  year    = {2026},
  url     = {https://github.com/<your-username>/biochemist-python},
  note    = {A complete workshop from Python basics to MD simulation}
}
```

---

## 📚 References

- [RCSB PDB](https://www.rcsb.org/) — Protein Data Bank
- [MolSSI](https://molssi.org/) — Molecular Sciences Software Institute
- [AutoDock Vina](https://vina.scripps.edu/) — Molecular Docking
- [ProLIF](https://prolif.readthedocs.io/) — Protein-Ligand Interaction Fingerprints
- [OpenMM](https://openmm.org/) — Molecular Dynamics Engine
- [RDKit](https://www.rdkit.org/) — Cheminformatics Toolkit

---

## 👤 Author

**Omar Hamouda**  
Computational Chemist | MolSSI Workshop Participant

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
