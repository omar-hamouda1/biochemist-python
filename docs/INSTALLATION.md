# ⚙️ Installation Guide

Complete setup instructions for running the Biochemist Python computational drug-discovery pipeline.

---

## Prerequisites

| Requirement | Minimum             | Recommended  |
| ----------- | ------------------- | ------------ |
| OS          | Ubuntu 20.04 / WSL2 | Ubuntu 22.04 |
| Python      | 3.10                | 3.11         |
| RAM         | 8 GB                | 16 GB        |
| Disk        | 5 GB                | 10 GB        |
| Conda       | Miniconda3          | Miniconda3   |

---

## Step 1: Install Miniconda

If Conda is not already installed:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

---

## Step 2: Clone the Repository

```bash
git clone https://github.com/omar-hamouda1/biochemist-python.git
cd biochemist-python
```

---

## Step 3: Create the Conda Environment

The recommended environment is defined in `environment.yml`:

```bash
conda env create -f environment.yml
conda activate biochem
```

This environment includes the main dependencies used by the project:

* **Core:** Python 3.11, NumPy, Pandas, SciPy, Matplotlib, Seaborn
* **Structure:** Biopython, MDAnalysis
* **Cheminformatics:** RDKit, Open Babel, RCSB search tools
* **Docking:** Smina, Meeko
* **Simulation:** OpenMM
* **Visualization:** NGLView, py3Dmol
* **Analysis:** ProLIF and related structural-analysis tools

For reproducible project setup, prefer `environment.yml` over installing individual packages manually.

---

## Step 4: Verify the Docking Environment

The current standardized docking workflow uses **Smina**.

Verify the executable:

```bash
smina --version
```

The standardized docking protocol is implemented in:

```text
scripts/run_standardized_docking.py
```

Its current protocol is:

* Receptor: `docking/receptor/2zq2_receptor.pdbqt`
* Center: `(17.672, -8.256, 10.688)`
* Box size: `25 × 25 × 25 Å`
* Exhaustiveness: `4`
* Number of modes: `1`
* Seed: `42`

The configuration is documented in:

```text
configs/docking_config.yml
```

The current standardized docking results are written to:

```text
docking/results/standardized/
docking/results/standardized_affinities.csv
```

The validated ligand manifest is:

```text
docking/results/validated_ligands_111.csv
```

---

## Step 5: Verify Open Babel

Open Babel is used for ligand-format conversion and structure preparation.

Verify the installation:

```bash
obabel -V
```

---

## Step 6: Launch JupyterLab

```bash
jupyter lab
```

The repository contains educational notebooks covering the project workflow.

The notebooks should be treated as educational and analysis material rather than as a required sequential execution path for reproducing the current standardized screening results.

The authoritative current docking workflow is the script:

```text
scripts/run_standardized_docking.py
```

---

## VS Code + WSL Setup

If you are using VS Code with WSL:

1. Install the **Remote - WSL** extension in VS Code.
2. Open VS Code and choose **Remote-WSL: Open Folder**.
3. Open the project directory, for example:

```text
/home/<user>/biochemist-python
```

4. Install the **Python** and **Jupyter** extensions.
5. Select the `biochem` Conda environment as the Python interpreter/kernel.

---

## Running Tests

Activate the project environment:

```bash
conda activate biochem
```

Run the complete test suite:

```bash
pytest -q
```

The repository currently contains 19 automated tests covering the core reusable modules and standardized docking utilities.

To check Python syntax for the project scripts:

```bash
python -m py_compile scripts/*.py
```

---

## Current Results and Reproducibility

The current standardized screening uses 117 candidate ligands as input.

After validation and resolution checks, 111 unique ligands are included in the standardized docking set.

Six ligands remain unresolved:

```text
0ZW
0ZX
0ZY
PPB
BAZ
BOZ
```

The authoritative current docking report is:

```text
docking/results/standardized_affinities.csv
```

The current Top 5 standardized docking results are:

```text
R11   -9.7876091
13U   -9.51114082
BAH   -9.41380882
12U   -9.3740406
T87   -9.27052689
```

These docking scores are computational ranking scores and should not be interpreted as experimental binding affinities.

Current ProLIF and rule-based drug-likeness results for the Top 5 are stored under:

```text
docking/results/prolif_summary.csv
docking/results/admet/
```

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'rdkit'`

Make sure the project environment is active:

```bash
conda activate biochem
```

Then verify:

```bash
python -c "from rdkit import Chem; print('RDKit OK')"
```

### `smina: command not found`

Make sure the `biochem` Conda environment is active:

```bash
conda activate biochem
```

Then verify:

```bash
smina --version
```

If Smina is still unavailable, recreate the environment from:

```bash
conda env create -f environment.yml
```

### `nglview` widget not showing

For notebooks where interactive widgets are unavailable, use the existing py3Dmol-based visualization cells when applicable.

### Open Babel conversion errors

Verify that Open Babel is available:

```bash
obabel -V
```

Also confirm that the input structure is valid and that the expected ligand/receptor files are present.

### MD simulation data too large

Large MD trajectory files are not committed to the repository.

Analysis notebooks may depend on locally available trajectory/analysis data under:

```text
md/colab_workshop/
```

Use the project documentation and existing analysis files to determine which local MD inputs are required for a specific notebook.
