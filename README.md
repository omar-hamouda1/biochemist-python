# 🧬 Biochemist Python — Computational Drug Discovery Pipeline

> **A hands-on computational drug discovery workshop for biochemists and life-science learners**
> From Python fundamentals and structural bioinformatics to virtual screening, ADMET analysis, and molecular dynamics.

**Target:** Trypsin (`2ZQ2`) — Serine Protease
**Reference ligand:** `13U`
**Screening library:** 117 ligands
**Validated docking set:** 111 ligands
**Documented exceptions:** 6 ligands

[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue?logo=python\&logoColor=white)](https://www.python.org/downloads/release/python-3110/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange?logo=jupyter)](https://jupyter.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2023.09.6-green)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-lightgrey.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-37%20passed-brightgreen)](tests/)
[![Conda](https://img.shields.io/badge/conda-environment-green?logo=anaconda)](environment.yml)

---

## 📌 Overview

**Biochemist Python** is a complete educational computational drug discovery pipeline built around a real protein target and a real screening workflow.

The project is designed to help biochemists learn how computational methods fit together in a practical research pipeline, including:

* Python programming and data analysis
* Protein and ligand structure handling
* Molecular visualization
* Ligand discovery from structural databases
* Binding-site analysis
* Standardized molecular docking
* Rule-based ADMET filtering
* Protein–ligand interaction analysis
* Molecular dynamics simulation and trajectory analysis

The project combines **educational Jupyter notebooks** with reusable Python modules, automated tests, configuration-driven docking, and dedicated audit utilities.

---

## 🎯 Current Scientific Workflow

| Item                      | Current project state        |
| ------------------------- | ---------------------------- |
| **Protein target**        | Trypsin — PDB `2ZQ2`         |
| **Reference ligand**      | `13U`                        |
| **Screening library**     | 117 ligands                  |
| **Validated docking set** | 111 ligands                  |
| **Documented exceptions** | 6 ligands                    |
| **Docking engine**        | Smina                        |
| **Docking configuration** | `configs/docking_config.yml` |
| **Docking audit**         | `scripts/audit_docking.py`   |
| **Automated tests**       | 37 passing                   |

The standardized screening workflow is configuration-driven and uses a single authoritative docking configuration rather than duplicating docking parameters across multiple scripts. Exact validation environments are captured separately by the provenance workflow; `environment.yml` is the reproducible environment specification, while the captured environment snapshot records the concrete packages used for a particular validated run.

---

## 🗺️ Pipeline

```text
Python Foundations
        │
        ▼
Data Analysis
        │
        ▼
Structural Bioinformatics
        │
        ▼
3D Visualization
        │
        ▼
Ligand Discovery
      (117)
        │
        ▼
Binding-Site Analysis
        │
        ▼
Standardized Virtual Screening
        │
        ├── Smina
        ├── 111 validated ligands
        ├── 6 documented exceptions
        └── automated audit
        │
        ▼
ADMET / Rule-Based Filtering
        │
        ▼
Protein–Ligand Interaction Analysis
        │
        ▼
Molecular Dynamics
        │
        ▼
Trajectory / Stability Analysis
```

For the detailed stage-by-stage workflow, see [`docs/PIPELINE.md`](docs/PIPELINE.md).

---

## ⚡ Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/omar-hamouda1/biochemist-python.git
cd biochemist-python
```

### 2. Create the Conda environment

```bash
conda env create -f environment.yml
```

Activate the environment defined by your local `environment.yml`.

The environment created by `environment.yml` is named `biochem`:

```bash
conda activate biochem
```

A local development environment may use a different name; the repository specification remains `biochem`.

### 3. Launch JupyterLab

```bash
jupyter lab
```

Then open:

```text
notebooks/
```

and begin with:

```text
01_python_basics.ipynb
```

For complete installation instructions, including VS Code, WSL, external tools, and environment setup, see [`docs/INSTALLATION.md`](docs/INSTALLATION.md).

---

## 📓 Educational Notebooks

| #  | Notebook                                                                               | Topic                                            | Stage         |
| -- | -------------------------------------------------------------------------------------- | ------------------------------------------------ | ------------- |
| 01 | [`01_python_basics.ipynb`](notebooks/01_python_basics.ipynb)                           | Python fundamentals for biochemists              | Foundations   |
| 02 | [`02_file_parsing.ipynb`](notebooks/02_file_parsing.ipynb)                             | Parsing CSV and PDB files                        | Foundations   |
| 03 | [`03_multiple_files.ipynb`](notebooks/03_multiple_files.ipynb)                         | Processing multiple files with `glob`            | Foundations   |
| 04 | [`04_pandas.ipynb`](notebooks/04_pandas.ipynb)                                         | Data analysis with Pandas                        | Foundations   |
| 05 | [`05_linear_regression.ipynb`](notebooks/05_linear_regression.ipynb)                   | Bradford protein assay — linear fit              | Data Analysis |
| 06 | [`06_plots.ipynb`](notebooks/06_plots.ipynb)                                           | Publication-quality plots                        | Data Analysis |
| 07 | [`07_nonlinear_regression_part1.ipynb`](notebooks/07_nonlinear_regression_part1.ipynb) | Michaelis–Menten kinetics                        | Data Analysis |
| 08 | [`08_nonlinear_regression_part2.ipynb`](notebooks/08_nonlinear_regression_part2.ipynb) | Enzyme inhibition kinetics                       | Data Analysis |
| 09 | [`09_biopython_mmcif.ipynb`](notebooks/09_biopython_mmcif.ipynb)                       | Structural bioinformatics with Biopython         | Structure     |
| 10 | [`10_rcsb_web_api.ipynb`](notebooks/10_rcsb_web_api.ipynb)                             | RCSB PDB REST API queries                        | Structure     |
| 11 | [`11_icn3d_visualization.ipynb`](notebooks/11_icn3d_visualization.ipynb)               | 3D structure visualization                       | Structure     |
| 12 | [`12_ec_class_ligands.ipynb`](notebooks/12_ec_class_ligands.ipynb)                     | Ligand discovery by EC class — 117 ligands       | Discovery     |
| 13 | [`13_binding_site.ipynb`](notebooks/13_binding_site.ipynb)                             | Binding-site analysis with MDAnalysis and ProLIF | Discovery     |
| 14 | [`14_molecular_docking.ipynb`](notebooks/14_molecular_docking.ipynb)                   | Molecular docking concepts and workflow          | Docking       |
| 15 | [`15_admet_prediction.ipynb`](notebooks/15_admet_prediction.ipynb)                     | ADMET / Lipinski Ro5 filtering                   | ADMET         |
| 16 | [`16_md_analysis.ipynb`](notebooks/16_md_analysis.ipynb)                               | MD simulation and trajectory analysis            | MD / Analysis |

> **Important:** The notebooks are educational material. The **authoritative standardized virtual-screening workflow** is the Python docking module and its validated configuration, not a notebook execution order.

---

## 🎯 Authoritative Standardized Docking Workflow

The production-standardized docking workflow is:

```text
scripts/run_standardized_docking.py
```

Its authoritative runtime configuration is:

```text
configs/docking_config.yml
```

The configuration defines the docking box, Smina settings, expected ligand count, input/output paths, manifest, and report location.

### Run standardized docking

From the repository root:

```bash
python -m scripts.run_standardized_docking
```

This workflow is resumable: valid existing docking outputs are reused rather than recomputed.

### Record reproducibility provenance

After preparing the receptor and ligands and completing standardized docking, capture the exact local environment and build the tracked provenance record:

```bash
python -m scripts.capture_environment
python -m scripts.build_provenance_manifest
python -m scripts.audit_docking --check-artifacts --check-provenance
```

The provenance manifest records configuration, receptor, prepared-ligand, docking-output, code, tool-version, and environment hashes. Generated PDBQT sidecars remain local integrity guards; the tracked manifest is the reviewable run-level record.

### Audit the standardized results

Report-level audit:

```bash
python -m scripts.audit_docking
```

Local artifact integrity audit:

```bash
python -m scripts.audit_docking --check-artifacts
```

The validated project state currently passes both audits:

```text
Candidate ligands       : 117
Validated ligands       : 111
Documented exceptions   : 6
Reported ligands        : 111
Expected validated     : 111
Status                  : PASS
```

The six documented exceptions are kept outside the validated screening set rather than being silently included or discarded.

---

## 🧪 Docking Configuration

The standardized protocol is defined in:

```text
configs/docking_config.yml
```

Current key parameters include:

```yaml
target_pdb: 2ZQ2
reference_ligand: 13U

center_x: 17.672
center_y: -8.256
center_z: 10.688

size_x: 25.0
size_y: 25.0
size_z: 25.0

smina_executable: smina
exhaustiveness: 4
num_modes: 1
seed: 42
timeout_seconds: 300

expected_ligands: 111
```

The former standalone box configuration file:

```text
docking/box_config.txt
```

has been retired in favor of the single validated YAML configuration.

---

## 🏆 Standardized Virtual Screening Results

The current standardized screening ranks the following compounds highest by docking score:

| Rank | Ligand  | Docking score (kcal/mol) | Lipinski   | Veber | PAINS |
| ---: | ------- | -----------------------: | ---------- | ----- | ----- |
|    1 | **R11** |                **-9.83** | Excellent  | Pass  | Clean |
|    2 | **13U** |                **-9.44** | Excellent  | Pass  | Clean |
|    3 | **BAH** |                **-9.42** | Acceptable | Fail  | Clean |
|    4 | **12U** |                **-9.34** | Excellent  | Pass  | Clean |
|    5 | **607** |                **-9.33** | Excellent  | Fail  | Clean |

These values are **docking scores**, not experimental binding measurements.

Docking results should therefore be interpreted as computational ranking evidence rather than direct measurements of biochemical potency or affinity.

---

## 🔬 Binding-Site and Interaction Analysis

Notebook 13 performs binding-site and protein–ligand interaction analysis for the Trypsin–`13U` reference system.

The validated structure-preparation workflow includes:

* deterministic protein alternate-location resolution via `scripts/prepare_protein.py`
* PDB2PQR-based protein protonation at pH 7.4 using PARSE
* RDKit ligand bond-order assignment
* Open Babel hydrogen addition
* ProLIF interaction analysis

The current validated ProLIF analysis identifies interactions involving residues including:

```text
ASN97.A
ASP189.A
GLN175.A
GLN192.A
GLY216.A
GLY219.A
LEU99.A
SER190.A
SER195.A
SER214.A
THR98.A
TRP215.A
```

This analysis is intended to support structural interpretation of the docking results rather than to replace experimental validation.

---

## 🧬 Molecular Dynamics Highlights

The MD results in this repository describe the historical **13U–Trypsin reference-ligand system**.

They should **not** be interpreted as molecular-dynamics validation of the current standardized docking top hit `R11`.

### System

* Complex: Trypsin–13U
* Solvent: explicit TIP3P water
* Salt: 0.15 M NaCl
* System size: 24,590 atoms
* Box size: 64.2 Å

### Simulation

* Production length: 5 ns
* Hardware: NVIDIA Tesla T4 GPU
* Reported performance: 196 ns/day

### Stability

* Protein backbone RMSD: **0.84 Å**
* Ligand RMSD: **1.57 Å**
* Ligand RMSF: **1.07 Å**
* Ligand atoms exceeding 3 Å RMSF: **0 / 65**

### MM/GBSA-style Analysis

The historical calculation uses **MMPBSA.py with the Generalized Born model (`igb=5`)** rather than a Poisson–Boltzmann model.

* ΔG estimate: **−32.13 ± 17.22 kcal/mol**
* van der Waals contribution: **−34.32 kcal/mol**
* Electrostatic contribution: **−11.22 kcal/mol**

This is a model-based computational estimate, not an experimental binding free energy.

### Key interaction residues

For the historical 13U–Trypsin system, ProLIF identified high-occupancy interactions involving:

* ASP171
* SER172
* SER192
* GLY196
* **TRP193** — largest favorable protein-side contribution in the stored decomposition

These MD and MM/GBSA-style results are historical/reference-system results and should be interpreted separately from the current 111-ligand standardized docking campaign. They do not experimentally validate the docking ranking.

---

## 🧪 Testing and Validation

The repository currently contains **37 automated tests**, all passing in the validated development environment:

```bash
pytest -q
```

Current result:

```text
37 passed
```

The test suite covers:

* docking output parsing
* standardized docking utilities
* validated YAML configuration loading
* path resolution
* Smina command construction
* docking-result auditing
* validated-manifest consistency
* documented exceptions
* invalid or incomplete standardized reports

The standardized docking audit also provides an independent project-level consistency check.

---

## 🛠️ Technology Stack

| Category                 | Tools                                        |
| ------------------------ | -------------------------------------------- |
| **Programming**          | Python 3.11                                  |
| **Data Analysis**        | NumPy, Pandas, SciPy                         |
| **Structure Analysis**   | Biopython, MDAnalysis                        |
| **Cheminformatics**      | RDKit, rcsb-api                         |
| **Docking**              | Smina, Open Babel             |
| **Interaction Analysis** | ProLIF                                       |
| **Visualization**        | py3Dmol, NGLView, iCN3D, Matplotlib, Seaborn |
| **Molecular Dynamics**   | OpenMM                                       |
| **ADMET / Filtering**    | RDKit descriptors and rule-based filters     |
| **Testing**              | pytest                                       |

---

## 📁 Project Structure

```text
biochemist-python/
│
├── notebooks/                    # Educational Jupyter notebooks
│   └── solutions/                # Solved workshop versions
│
├── src/                          # Reusable Python modules
│   ├── protein.py                # Protein structure utilities
│   ├── admet.py                  # ADMET / Lipinski calculations
│   ├── docking.py                # Docking result parsing/utilities
│   ├── docking_config.py         # Validated docking configuration model
│   ├── provenance.py             # Artifact provenance and integrity helpers
│   └── visualization.py          # 3D visualization helpers
│
├── scripts/                      # Standalone workflow scripts
│   ├── audit_docking.py          # Standardized docking audit
│   ├── run_standardized_docking.py
│   ├── prepare_ligands.py        # Provenance-aware ligand preparation
│   ├── prepare_protein.py        # Deterministic protein altLoc/PDB2PQR preparation
│   ├── prepare_receptor.py       # Deterministic receptor PDBQT preparation
│   ├── build_provenance_manifest.py
│   └── capture_environment.py    # Exact validation-environment snapshot
│
├── tests/                        # Automated tests
│
├── configs/                      # Authoritative configuration files
│   └── docking_config.yml
│
├── ligands/                      # 117 ligand SDF files
│
├── docking/                      # Docking inputs, manifests, results
│   ├── ligands/
│   ├── receptor/
│   └── results/
│
├── pdb/                          # Protein and ligand structure files
├── pdb_files/                    # Additional structural datasets
├── md/                           # Molecular dynamics data
├── figures/                      # Generated figures
├── results/                      # Analysis results
├── data/                         # Input datasets
├── molssi_data/                  # MolSSI workshop material
│
├── docs/                         # Project documentation
│   ├── PIPELINE.md
│   └── INSTALLATION.md
│
├── environment.yml               # Conda environment specification
├── requirements.txt              # Python package requirements
├── LICENSE                       # MIT License
├── CONTRIBUTING.md               # Contribution guidelines
├── CODE_OF_CONDUCT.md            # Community standards
└── CHANGELOG.md                  # Project changelog
```

---

## 💾 Data Availability

Large molecular-dynamics trajectory files, including:

```text
*.dcd
*.nc
```

are not included in the repository when they exceed practical GitHub size limits.

For the full MD analysis workflow, download the required trajectory data and place it in:

```text
md/colab_workshop/02_analysis/
```

---

## 📚 Documentation

Detailed documentation is available in:

* [`docs/INSTALLATION.md`](docs/INSTALLATION.md) — installation and environment setup
* [`docs/PIPELINE.md`](docs/PIPELINE.md) — computational workflow and methodology

---

## 🤝 Contributing

Contributions, corrections, documentation improvements, and scientific feedback are welcome.

Please read [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening a pull request.

---

## 📚 Citation

If you use this project in teaching, research, or derived computational work, please cite:

```bibtex
@software{hamouda2025biochemist,
  author  = {Hamouda, Omar},
  title   = {Biochemist Python: Computational Drug Discovery Pipeline},
  year    = {2025},
  url     = {https://github.com/omar-hamouda1/biochemist-python},
  note    = {A hands-on computational drug discovery workshop from Python fundamentals to molecular dynamics}
}
```

---

## 📖 References

* [RCSB Protein Data Bank](https://www.rcsb.org/)
* [AutoDock Vina](https://vina.scripps.edu/)
* [ProLIF](https://prolif.readthedocs.io/)
* [OpenMM](https://openmm.org/)
* [RDKit](https://www.rdkit.org/)

---

## 👤 Author

**Omar Hamouda**
B.Sc. Chemistry, Faculty of Science — Suez Canal University

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
