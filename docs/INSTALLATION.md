# ⚙️ Installation Guide

Complete setup instructions for running the Biochemist Python pipeline.

---

## Prerequisites

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| OS | Ubuntu 20.04 / WSL2 | Ubuntu 22.04 |
| Python | 3.10 | 3.11 |
| RAM | 8 GB | 16 GB |
| Disk | 5 GB | 10 GB |
| Conda | Miniconda3 | Miniconda3 |

---

## Step 1: Install Miniconda (if not installed)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

## Step 2: Clone the Repository

```bash
git clone https://github.com/omar-hamouda1/biochemist-python.git
cd biochemist-python
```

## Step 3: Create the Conda Environment

```bash
conda env create -f environment.yml
conda activate biochem
```

This installs all dependencies:
- **Core**: Python 3.11, NumPy, Pandas, SciPy, Matplotlib, Seaborn
- **Structure**: Biopython, MDAnalysis
- **Cheminformatics**: RDKit, rcsbsearchapi
- **Docking**: (requires separate AutoDock Vina installation)
- **Simulation**: OpenMM
- **Visualization**: NGLView, ICN3D, py3Dmol
- **Analysis**: ProLIF, PDB2PQR

## Step 4: Install AutoDock Vina (for Notebook 14)

```bash
# Option A: via conda
conda install -c conda-forge vina

# Option B: via apt (Ubuntu)
sudo apt install autodock-vina

# Option C: Download binary
wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_2_5_linux_x86.tgz
tar xzf autodock_vina_1_2_5_linux_x86.tgz
sudo cp bin/vina /usr/local/bin/
```

Verify:
```bash
vina --version
```

## Step 5: Install Open Babel (for Notebook 14)

```bash
conda install -c conda-forge openbabel
```

## Step 6: Launch JupyterLab

```bash
jupyter lab
```

Navigate to `notebooks/` and run notebooks **01 → 16** in order.

---

## VS Code + WSL Setup

If you're using VS Code with WSL:

1. Install the **Remote - WSL** extension in VS Code
2. Open VS Code → `Ctrl+Shift+P` → "Remote-WSL: Open Folder"
3. Navigate to `/home/<user>/biochemist-python`
4. Install the **Jupyter** extension in VS Code
5. Select the `biochem` conda environment as kernel

---

## Running Tests

```bash
conda activate biochem
pytest tests/ -v
```

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'rdkit'"
→ Make sure you activated the conda env: `conda activate biochem`

### "nglview widget not showing"
→ In VS Code, use py3Dmol instead (already configured in notebooks)

### "vina: command not found"
→ Install AutoDock Vina (see Step 4 above)

### MD Simulation data too large
→ The `.dcd` trajectory files are not included in the repo.
   Download them from [Zenodo link] and place in `md/colab_workshop/02_analysis/`
