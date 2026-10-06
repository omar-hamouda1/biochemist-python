# Docking Results

## Authoritative current results

These files represent the current standardized virtual-screening pipeline:

- `standardized_affinities.csv` — authoritative 111-ligand docking affinity report.
- `validated_ligands_111.csv` — manifest of the 111 validated ligand IDs used by the standardized docking script.
- `admet/admet_summary.csv` — current rule-based descriptor and filter summary for the current Top 5.
- `admet/admet_top5.csv` — current Top 5 molecular-property/drug-likeness summary.
- `prolif_summary.csv` — current ProLIF summary for the current Top 5.
- `ligand_preparation.csv` — ligand preparation status for the input set.

The current standardized docking protocol is defined in:

`configs/docking_config.yml`

and executed by:

`scripts/run_standardized_docking.py`

The standardized docking protocol uses Smina with:

- receptor: `docking/receptor/2zq2_receptor.pdbqt`
- center: `(17.672, -8.256, 10.688)`
- box size: `25 × 25 × 25 Å`
- exhaustiveness: `4`
- num modes: `1`
- seed: `42`

## Recovery / intermediate results

These files belong to earlier recovery or intermediate stages and are not the authoritative current ranking:

- `affinities_recovered_111.csv`
- `recovered_affinities.csv`
- `docking_exceptions.csv`

`affinities_recovered_111.csv` was used to establish the validated 111-ligand set during the recovery stage. Its affinity values are not the current standardized docking scores.

## Historical / educational outputs

The following files are retained for historical or educational context and should not be used as the current docking ranking:

- `affinities_full.csv`
- `affinities.csv`
- `prolif_final_summary.csv`
- `top_5.txt`
- older ADMET component files such as `admet_descriptors.csv`, `admet_lipinski.csv`, `admet_veber.csv`, and `admet_pains.csv`
- older fixed-structure, complex, and visualization outputs retained from earlier analysis stages

When current results are needed, use the authoritative files listed above.