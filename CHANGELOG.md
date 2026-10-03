# 📋 Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.0.0] - 2026-10-02

### Added
- Complete 16-notebook pipeline (Python basics → MD Simulation)
- Virtual screening of 117 ligands against Trypsin (2ZQ2)
- AutoDock Vina docking with automated batch processing
- ADMET prediction module () with Lipinski Ro5
- ProLIF interaction fingerprint analysis ()
- MD simulation setup with OpenMM (minimization, equilibration, production)
- Binding site investigation with py3Dmol / NGLView
- Unit tests for ADMET module
- Comprehensive documentation (README, CONTRIBUTING, CODE_OF_CONDUCT)
- GitHub Actions CI pipeline
- Professional project structure with , , 
- Environment configuration () with pinned versions
- Docking box configuration for reproducibility

### Data
- 117 ligand SDF files from PDB (EC class 3.4.21 — Serine Proteases)
- Docking results:  with ranked binding energies
- Top hits: 13U (-9.52), R11 (-9.51), BAH (-9.40), 607 (-9.40), 12U (-9.35)
- ProLIF interaction analysis for top 4 ligands
- MD simulation figures: RMSD, RMSF, H-bonds, MMPBSA decomposition

## [0.1.0] - 2026-08-15

### Added
- Initial project setup from MolSSI workshop materials
- Basic notebook structure (01-08: Python fundamentals)
- Protein structure files (2ZQ2.pdb)
