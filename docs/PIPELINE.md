# 🔬 Pipeline Architecture — Detailed Walkthrough

This document explains the computational drug-discovery workflow implemented in this project. Each stage maps to one or more Jupyter notebooks, reusable modules, and/or standalone scripts.

---

## Pipeline Overview

```text
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
Stage 6: Analysis            Stage 5: Properties         Stage 4: Docking
┌──────────────────┐         ┌──────────────────┐       ┌──────────────────┐
│ MD Trajectory    │         │ Lipinski / Veber │       │ Standardized     │
│ RMSD / RMSF      │◀────────│ Drug-likeness    │◀──────│ Smina            │
│ MM/GBSA-style   │         │ PAINS Filtering  │       │ Virtual Screen   │
│ ProLIF Analysis  │         │ RDKit Descriptors│       │ 111 ligands      │
└──────────────────┘         └──────────────────┘       │ Ranked by score  │
 NB: 16                       NB: 15                    └──────────────────┘
                                                               NB: 14
```

---

## Stage 1: Python Foundations (Notebooks 01–08)

### What You Learn

* Python syntax, variables, loops, functions, and conditionals
* File I/O: reading CSV and PDB files
* Processing multiple files with `glob` and `os`
* Data manipulation with `pandas`
* Linear regression (Bradford protein assay)
* Plotting with `matplotlib` and `seaborn`
* Nonlinear regression (Michaelis-Menten enzyme kinetics)

### Key Concepts

| Notebook | Biochemistry Concept | Python Concept             |
| -------- | -------------------- | -------------------------- |
| 01       | —                    | Variables, types, loops    |
| 02       | PDB file format      | String parsing, file I/O   |
| 03       | Batch analysis       | `glob`, `os.path`          |
| 04       | Data tables          | `pandas` DataFrames        |
| 05       | Bradford assay       | `scipy.stats.linregress`   |
| 06       | Assay visualization  | `matplotlib`, `seaborn`    |
| 07       | Michaelis-Menten     | `scipy.optimize.curve_fit` |
| 08       | Inhibition kinetics  | Nonlinear regression       |

---

## Stage 2: Structural Bioinformatics (Notebooks 09–11)

### What You Learn

* Parsing mmCIF files with Biopython
* Querying the RCSB PDB
* 3D structure visualization with ICN3D and py3Dmol
* Understanding resolution, R-factor, and structure quality

### Key Data

* **Myoglobin structures:** 40+ CIF files in `pdb_files/`
* **Target protein:** Trypsin (2ZQ2) — 1.40 Å resolution

---

## Stage 3: Drug Discovery (Notebooks 12–13)

### What You Learn

* Searching PDB by EC classification (EC 3.4.21 → Serine Proteases)
* Extracting bound ligands from crystal structures
* Building a screening library of 117 ligand SDF files
* Binding-site identification and visualization
* Understanding protein-ligand interactions

### Output

* `ligands/` directory: 117 ligand SDF files
* Standardized screening set: 111 validated ligands
* Binding-site residue information
* Reference ligand (13U) interaction profile

---

## Stage 4: Standardized Virtual Screening

### What You Learn

* Receptor preparation (PDB → PDBQT)
* Ligand preparation (SDF → PDBQT)
* Docking-box definition
* Running standardized Smina docking
* Parsing and ranking docking results
* Distinguishing current authoritative results from historical recovery outputs

### Authoritative Workflow

The current standardized screening is implemented by:

```text
scripts/run_standardized_docking.py
```

The authoritative runtime configuration is:

```text
configs/docking_config.yml
```

The docking script loads and validates this YAML at runtime; production docking parameters and paths are not duplicated as independent constants in the script.

Run the authoritative workflow from the repository root as:

```bash
python -m scripts.run_standardized_docking
```

The companion audit for the authoritative report is:

```bash
python -m scripts.audit_docking
```

Use `--check-artifacts` with the audit when generated receptor, ligand, and docking PDBQT files are available locally.

The validated ligand manifest is:

```text
docking/results/validated_ligands_111.csv
```

### Configuration

```text
Receptor:  docking/receptor/2zq2_receptor.pdbqt
Center:    (17.672, -8.256, 10.688) Å
Box size:  25.0 × 25.0 × 25.0 Å
Tool:      Smina
Exhaustiveness: 4
Number of modes: 1
Random seed: 42
```

The standardized run is defined for the validated 111-ligand set. The expected ligand count, manifest, report path, timeout, and Smina executable are controlled by `configs/docking_config.yml`.

The authoritative docking report is:

```text
docking/results/standardized_affinities.csv
```

### Receptor and ligand provenance

Protein preparation is deterministic and recorded by:

```bash
python -m scripts.prepare_protein
```

This resolves protein alternate locations by highest mean occupancy per residue, writes a coordinate-only PDB for PDB2PQR, and runs PDB2PQR at pH 7.4 with PARSE.

The prepared PDB2PQR product is then converted to the rigid receptor PDBQT by:

```bash
python -m scripts.prepare_receptor
```

It converts the PDB2PQR-prepared `pdb/protein_h.pdb` to the rigid receptor PDBQT using the documented Open Babel protocol.

Ligand preparation is project-rooted and provenance-aware:

```bash
python -m scripts.prepare_ligands
```

Each generated PDBQT receives a local provenance sidecar. These sidecars are intentionally ignored by Git; the reviewable run-level record is generated after the full screening run:

```bash
python -m scripts.capture_environment
python -m scripts.build_provenance_manifest
python -m scripts.audit_docking --check-artifacts --check-provenance
```

The tracked `docking/results/standardized_provenance.json` records hashes for the configuration, source structures, protein/receptor preparation and docking code, tools, prepared ligands, docking poses, and validation environment.

The standardized pose files are written under:

```text
docking/results/standardized/
```

Historical recovery-stage affinity files are retained separately and are not used as the current docking ranking.

### Current Top 5 Results

| Rank | Ligand | Affinity (kcal/mol) |
| ---- | ------ | ------------------- |
| 1    | R11    | -9.834149           |
| 2    | 13U    | -9.438036           |
| 3    | BAH    | -9.420844           |
| 4    | 12U    | -9.339968           |
| 5    | 607    | -9.332170           |

These are computational docking scores used for ranking and should not be interpreted as experimental binding affinities.

---

## Stage 5: Rule-Based Drug-Likeness and Molecular Properties (Notebook 15 + ADMET script)

### What You Learn

* Lipinski's Rule of Five
* Molecular descriptors such as MW, LogP, HBD, HBA, TPSA, and RotB
* Veber-rule assessment
* PAINS filtering
* Rule-based drug-likeness assessment
* RDKit descriptor-based screening

The current workflow uses RDKit descriptors and rule-based filters. It is not a full ADME or toxicity prediction system.

The current Top 5 summary is stored under:

```text
docking/results/admet/
```

The authoritative current files include:

```text
docking/results/admet/admet_summary.csv
docking/results/admet/admet_top5.csv
```

Notebook 15 is retained as educational analysis material and should not be treated as a separate authoritative ranking pipeline.

---

## Stage 6: Molecular Dynamics (Notebook 16)

### What You Learn

* System preparation (protein + ligand + solvent + ions)
* Energy minimization
* NVT and NPT equilibration
* Production MD
* Trajectory analysis (RMSD, RMSF, hydrogen bonds)
* MM/GBSA-style binding-energy decomposition using the Generalized Born model
* ProLIF interaction analysis

The MD analysis documented in this repository is based on the historical Trypsin–13U reference-ligand system.

The MD results provide model-based computational evidence about the behavior of 13U in that simulated complex. They do not experimentally validate the current docking ranking, and they do not constitute MD validation of the current top-ranked ligand R11.

### Software

* **OpenMM** for simulation
* **MDAnalysis** for trajectory analysis
* **ProLIF** for interaction fingerprints

---

## Current Results and Reproducibility

### Validated Screening Set

The project starts from 117 candidate ligands.

After validation and resolution checks:

* 111 unique ligands are included in the standardized docking set.
* 6 ligands remain unresolved: `0ZW`, `0ZX`, `0ZY`, `PPB`, `BAZ`, and `BOZ`.

### Authoritative Current Outputs

| Purpose                      | File                                          |
| ---------------------------- | --------------------------------------------- |
| Validated ligand manifest    | `docking/results/validated_ligands_111.csv`   |
| Standardized docking report  | `docking/results/standardized_affinities.csv` |
| ProLIF Top 5 summary         | `docking/results/prolif_summary.csv`          |
| Molecular-property Top 5 summary | `docking/results/admet/admet_top5.csv`        |
| Merged Top-hit summary       | `results/top_hits_summary.csv`                |

For the current results, prefer the authoritative files listed above.

---

## Directory Map

```text
biochemist-python_ORGANIZED/
├── notebooks/          # 16 educational/analysis notebooks
├── src/                # Reusable Python modules
│   ├── protein.py      # Protein structure utilities
│   ├── admet.py        # Drug-likeness/property calculations
│   ├── docking.py      # Standardized docking result parsing
│   ├── docking_config.py# Validated docking configuration model
│   └── visualization.py# 3D visualization helpers
├── scripts/            # Standalone analysis and pipeline scripts
├── tests/              # Unit tests (pytest)
├── data/               # Educational/reference data
├── pdb/                # Protein structure files
├── pdb_files/          # Myoglobin CIF collection
├── ligands/            # 117 ligand SDF files
├── docking/            # Docking inputs and results
├── md/                 # MD simulation files and analysis
├── figures/            # Generated figures
├── results/            # Final merged result summaries
├── configs/            # Configuration files
├── molssi_data/        # MolSSI workshop reference material
└── docs/               # Project documentation
```
