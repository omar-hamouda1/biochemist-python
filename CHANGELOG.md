# 📋 Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Changed

* Aligned the documented docking configuration with the standardized Smina protocol.
* Standardized docking utilities and tests to use `REMARK minimizedAffinity` outputs.
* Separated the validated 111-ligand manifest from historical docking recovery results.
* Clarified the distinction between current authoritative results, recovery/intermediate artifacts, and historical/educational outputs.
* Updated MD analysis wording to distinguish historical 13U simulation evidence from the current standardized docking ranking.
* Updated the docking results documentation to identify authoritative current outputs.

### Fixed

* Corrected legacy references to AutoDock Vina in the active docking utilities.
* Corrected outdated Top 5 and ProLIF references in the MD analysis notebook.
* Removed wording that implied experimental validation of the current docking ranking by the historical 13U MD simulation.

## [1.0.0] - 2026-10-02

### Added

* Complete 16-notebook educational pipeline covering Python fundamentals through molecular dynamics analysis.
* Virtual screening workflow for 117 candidate ligands against Trypsin (PDB: 2ZQ2).
* Standardized Smina docking workflow for the validated 111-ligand set.
* Rule-based drug-likeness and molecular-property assessment using Lipinski, Veber, and PAINS filters.
* ProLIF protein-ligand interaction fingerprint analysis.
* MD simulation workflow using OpenMM for minimization, equilibration, and production.
* Binding-site investigation with py3Dmol / NGLView.
* Unit tests and GitHub Actions CI pipeline.
* Comprehensive project documentation and reproducible environment configuration.
* Docking box configuration for reproducible standardized screening.

### Data

* 117 ligand SDF files from PDB-derived serine protease ligands.
* 111 validated unique ligands in the standardized docking set.
* 6 unresolved ligands excluded from the standardized docking set: `0ZW`, `0ZX`, `0ZY`, `PPB`, `BAZ`, and `BOZ`.
* Authoritative standardized docking report: `docking/results/standardized_affinities.csv`.
* Current Top 5 standardized docking results:

  * R11: `-9.7876091`
  * 13U: `-9.51114082`
  * BAH: `-9.41380882`
  * 12U: `-9.3740406`
  * T87: `-9.27052689`
* Current ProLIF interaction analysis for the standardized Top 5.
* Current rule-based drug-likeness and molecular-property summary for the standardized Top 5.
* Historical MD analysis of the Trypsin–13U complex, including RMSD, RMSF, hydrogen-bond, ProLIF, and MM/GBSA-style analysis using MMPBSA.py with the Generalized Born model.

## [0.1.0] - 2026-08-15

### Added

* Initial project setup from MolSSI workshop materials.
* Basic notebook structure (01–08: Python fundamentals).
* Protein structure files and initial Trypsin structure analysis (2ZQ2).
