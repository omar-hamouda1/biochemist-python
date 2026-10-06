"""Audit the authoritative standardized docking results."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Iterable

from src.docking import parse_smina_output
from src.docking_config import DockingConfig, load_docking_config
from src.provenance import (
    build_metadata,
    executable_version,
    load_metadata,
    provenance_path,
    sha256_file,
    validate_existing_output,
)


@dataclass(frozen=True)
class DockingAuditReport:
    """Structured audit results for the standardized docking pipeline."""

    candidate_ids: tuple[str, ...]
    validated_ids: tuple[str, ...]
    exception_ids: tuple[str, ...]
    report_ids: tuple[str, ...]
    errors: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.errors


def _sorted_unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted(set(values)))


def _read_csv_rows(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _read_id_column(path: Path, label: str) -> tuple[tuple[str, ...], list[str]]:
    rows = _read_csv_rows(path, label)
    if not rows or "ligand_id" not in rows[0]:
        raise ValueError(f"{label} must contain a 'ligand_id' column: {path}")
    raw_ids = [row.get("ligand_id", "").strip() for row in rows]
    invalid = [value for value in raw_ids if not value]
    if invalid:
        raise ValueError(f"{label} contains an empty ligand_id: {path}")
    duplicates = sorted({value for value in raw_ids if raw_ids.count(value) > 1})
    return _sorted_unique(raw_ids), duplicates


def _audit_report_csv(
    path: Path,
    expected_ids: set[str],
    errors: list[str],
) -> tuple[str, ...]:
    rows = _read_csv_rows(path, "standardized docking report")
    required = {"ligand_id", "affinity", "status"}
    if not rows:
        errors.append(f"Standardized docking report is empty: {path}")
        return ()

    missing_columns = required - set(rows[0])
    if missing_columns:
        errors.append(
            "Standardized docking report is missing columns: "
            + ", ".join(sorted(missing_columns))
        )
        return ()

    report_ids_raw = [row.get("ligand_id", "").strip() for row in rows]
    duplicates = sorted(
        {value for value in report_ids_raw if value and report_ids_raw.count(value) > 1}
    )
    if duplicates:
        errors.append(f"Duplicate ligand IDs in standardized report: {duplicates}")

    for row in rows:
        ligand_id = row.get("ligand_id", "").strip()
        status = row.get("status", "").strip()
        affinity_text = row.get("affinity", "").strip()

        if not ligand_id:
            errors.append("Standardized report contains an empty ligand_id")
            continue
        if status not in {"ok", "existing"}:
            errors.append(
                f"Standardized report has non-complete status for {ligand_id}: {status!r}"
            )
        try:
            affinity = float(affinity_text)
        except (TypeError, ValueError):
            errors.append(f"Non-numeric affinity for {ligand_id}: {affinity_text!r}")
            continue
        if not isfinite(affinity):
            errors.append(f"Non-finite affinity for {ligand_id}: {affinity_text!r}")

    report_ids = set(report_ids_raw)
    missing = sorted(expected_ids - report_ids)
    extra = sorted(report_ids - expected_ids)
    if missing:
        errors.append(f"Validated ligands missing from standardized report: {missing}")
    if extra:
        errors.append(f"Unknown ligands in standardized report: {extra}")

    return _sorted_unique(report_ids_raw)


def _audit_provenance(
    config: DockingConfig,
    validated_ids: Iterable[str],
    errors: list[str],
) -> None:
    """Audit provenance sidecars and the tracked run manifest."""
    manifest_path = config.report.parent / "standardized_provenance.json"
    if not manifest_path.exists():
        errors.append(f"Missing tracked provenance manifest: {manifest_path}")
        return

    try:
        import json

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"Invalid provenance manifest: {exc}")
        return

    if manifest.get("schema_version") != 1:
        errors.append("Unsupported provenance manifest schema")

    expected_config_sha = sha256_file(
        config.report.parents[2] / "configs" / "docking_config.yml"
    )
    if manifest.get("config", {}).get("sha256") != expected_config_sha:
        errors.append("Provenance manifest config hash does not match current config")

    receptor = config.receptor
    receptor_meta = load_metadata(provenance_path(receptor))
    if receptor_meta is None:
        errors.append(f"Missing receptor provenance: {provenance_path(receptor)}")
    elif receptor_meta.get("output_sha256") != sha256_file(receptor):
        errors.append("Receptor provenance output hash does not match receptor")

    runner_path = config.report.parents[2] / "scripts" / "run_standardized_docking.py"
    parser_path = config.report.parents[2] / "src" / "docking.py"
    smina_version = executable_version(config.smina_executable)

    for ligand_id in validated_ids:
        ligand = config.ligand_dir / f"{ligand_id}.pdbqt"
        output = config.output_dir / f"{ligand_id}_out.pdbqt"
        ligand_meta = load_metadata(provenance_path(ligand))

        if ligand_meta is None:
            errors.append(f"{ligand_id}: missing ligand provenance")
        else:
            if ligand_meta.get("source_sha256") != sha256_file(
                config.report.parents[2] / "ligands" / f"{ligand_id}.sdf"
            ):
                errors.append(f"{ligand_id}: ligand source hash mismatch")
            if not ligand.exists():
                errors.append(f"{ligand_id}: prepared ligand missing")
            elif ligand_meta.get("output_sha256") != sha256_file(ligand):
                errors.append(f"{ligand_id}: prepared ligand hash mismatch")

        output_meta = load_metadata(provenance_path(output))
        if output_meta is None:
            errors.append(f"{ligand_id}: missing docking provenance")
            continue
        if not output.exists():
            errors.append(f"{ligand_id}: docking output missing")
            continue

        affinities = parse_smina_output(str(output))
        if len(affinities) != 1:
            errors.append(
                f"{ligand_id}: expected exactly one docking affinity, found {len(affinities)}"
            )

        expected = build_metadata(
            config_path=config.report.parents[2] / "configs" / "docking_config.yml",
            receptor_path=receptor,
            ligand_path=ligand,
            runner_path=runner_path,
            parser_path=parser_path,
            smina_version=smina_version,
            command=[
                config.smina_executable,
                "--receptor", str(config.receptor),
                "--ligand", str(ligand),
                "--center_x", str(config.center_x),
                "--center_y", str(config.center_y),
                "--center_z", str(config.center_z),
                "--size_x", str(config.size_x),
                "--size_y", str(config.size_y),
                "--size_z", str(config.size_z),
                "--exhaustiveness", str(config.exhaustiveness),
                "--num_modes", str(config.num_modes),
                "--seed", str(config.seed),
                "--out", str(output),
            ],
        )
        if not validate_existing_output(
            output_path=output,
            metadata_path=provenance_path(output),
            expected_metadata=expected,
            expected_affinity=affinities[0] if len(affinities) == 1 else None,
        ):
            errors.append(f"{ligand_id}: docking provenance validation failed")

    try:
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        current_report_sha = sha256_file(config.report)
        recorded_report_sha = (
            manifest_data.get("authoritative_files", {})
            .get("standardized_report", {})
            .get("sha256")
        )
        if recorded_report_sha != current_report_sha:
            errors.append("Provenance manifest standardized-report hash is stale")
    except (OSError, json.JSONDecodeError):
        pass


def audit_standardized_docking(
    config: DockingConfig,
    *,
    exceptions_path: Path | None = None,
    check_artifacts: bool = False,
    check_provenance: bool = False,
) -> DockingAuditReport:
    """Audit the authoritative screening manifest and standardized report.

    The audit deliberately ignores historical/recovery CSVs. Generated PDBQT
    artifacts are optional unless ``check_artifacts=True`` because those files
    are intentionally ignored by Git and are not required for report-level
    consistency checks.
    """

    errors: list[str] = []
    warnings: list[str] = []

    project_root = config.manifest.parents[2]
    candidate_dir = project_root / "ligands"
    candidate_set_available = candidate_dir.exists()
    if candidate_set_available:
        candidate_ids = _sorted_unique(
            path.stem for path in candidate_dir.glob("*.sdf")
        )
    else:
        candidate_ids = ()
        warnings.append(
            "Raw ligand directory not found; candidate-set partition was not checked: "
            f"{candidate_dir}"
        )

    try:
        validated_ids, manifest_duplicates = _read_id_column(
            config.manifest, "validated ligand manifest"
        )
    except (FileNotFoundError, ValueError) as exc:
        errors.append(str(exc))
        validated_ids = ()
        manifest_duplicates = []

    if manifest_duplicates:
        errors.append(f"Duplicate IDs in validated manifest: {manifest_duplicates}")

    if len(validated_ids) != config.expected_ligands:
        errors.append(
            f"Expected {config.expected_ligands} validated ligands, "
            f"found {len(validated_ids)}"
        )

    candidate_set = set(candidate_ids)
    validated_set = set(validated_ids)

    if candidate_set_available:
        unknown_validated = sorted(validated_set - candidate_set)
        if unknown_validated:
            errors.append(
                f"Validated manifest contains unknown ligands: {unknown_validated}"
            )

    if exceptions_path is None:
        exceptions_path = config.report.parent / "docking_exceptions.csv"

    try:
        exception_rows = _read_csv_rows(exceptions_path, "docking exceptions")
        if not exception_rows or "ligand_id" not in exception_rows[0]:
            raise ValueError(
                f"Docking exceptions must contain a 'ligand_id' column: {exceptions_path}"
            )
        exception_ids_raw = [row.get("ligand_id", "").strip() for row in exception_rows]
        if any(not value for value in exception_ids_raw):
            raise ValueError(
                f"Docking exceptions contains an empty ligand_id: {exceptions_path}"
            )
        exception_ids = _sorted_unique(exception_ids_raw)
        duplicate_exceptions = sorted(
            {value for value in exception_ids_raw if exception_ids_raw.count(value) > 1}
        )
        if duplicate_exceptions:
            errors.append(f"Duplicate IDs in docking exceptions: {duplicate_exceptions}")
        for row in exception_rows:
            if row.get("status", "").strip() != "unresolved":
                errors.append(
                    "Docking exceptions contains a non-unresolved row for "
                    f"{row.get('ligand_id', '').strip()}: {row.get('status', '').strip()!r}"
                )
    except (FileNotFoundError, ValueError) as exc:
        errors.append(str(exc))
        exception_ids = ()

    exception_set = set(exception_ids)
    if candidate_set_available:
        unknown_exceptions = sorted(exception_set - candidate_set)
        if unknown_exceptions:
            errors.append(
                f"Docking exceptions contain unknown ligands: {unknown_exceptions}"
            )

    overlap = sorted(validated_set & exception_set)
    if overlap:
        errors.append(
            f"Ligands appear in both validated manifest and exceptions: {overlap}"
        )

    if candidate_set_available:
        expected_partition = validated_set | exception_set
        missing_partition = sorted(candidate_set - expected_partition)
        unexpected_partition = sorted(expected_partition - candidate_set)
        if missing_partition:
            errors.append(
                "Candidate ligands are neither validated nor documented as exceptions: "
                f"{missing_partition}"
            )
        if unexpected_partition:
            errors.append(
                "Manifest/exception IDs missing from ligand directory: "
                f"{unexpected_partition}"
            )

    report_ids = _audit_report_csv(
        config.report,
        validated_set,
        errors,
    )

    if check_artifacts:
        if not config.receptor.exists() or config.receptor.stat().st_size == 0:
            errors.append(f"Missing or empty receptor PDBQT: {config.receptor}")
        for ligand_id in validated_ids:
            ligand = config.ligand_dir / f"{ligand_id}.pdbqt"
            output = config.output_dir / f"{ligand_id}_out.pdbqt"
            if not ligand.exists() or ligand.stat().st_size == 0:
                errors.append(f"Missing or empty prepared ligand PDBQT: {ligand}")
            if not output.exists() or output.stat().st_size == 0:
                errors.append(f"Missing or empty docking output PDBQT: {output}")

    if check_provenance:
        _audit_provenance(config, validated_ids, errors)

    if candidate_set_available and len(candidate_ids) != len(validated_ids) + len(exception_ids):
        errors.append(
            "Candidate partition count mismatch: "
            f"{len(candidate_ids)} candidates != "
            f"{len(validated_ids)} validated + {len(exception_ids)} exceptions"
        )

    if not check_artifacts:
        warnings.append(
            "Generated receptor/ligand/docking PDBQT files were not checked; "
            "use --check-artifacts for a local artifact integrity audit."
        )

    return DockingAuditReport(
        candidate_ids=candidate_ids,
        validated_ids=validated_ids,
        exception_ids=exception_ids,
        report_ids=report_ids,
        errors=tuple(errors),
        warnings=tuple(warnings),
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit the authoritative standardized docking results."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to docking_config.yml (default: project config)",
    )

    parser.add_argument(
        "--exceptions",
        type=Path,
        default=None,
        help="Path to docking_exceptions.csv (default: beside the report)",
    )

    parser.add_argument(
        "--check-artifacts",
        action="store_true",
        help="Also check generated receptor, ligand, and docking PDBQT files.",
    )

    parser.add_argument(
        "--check-provenance",
        action="store_true",
        help="Also validate provenance manifests and docking sidecars.",
    )

    args = parser.parse_args()
    config = load_docking_config(args.config) if args.config else load_docking_config()
    audit = audit_standardized_docking(
        config,
        exceptions_path=args.exceptions,
        check_artifacts=args.check_artifacts,
        check_provenance=args.check_provenance,
    )

    print("=" * 80)
    print("STANDARDIZED DOCKING AUDIT")
    print("=" * 80)
    print(f"Target             : {config.target_pdb}")
    print(f"Reference ligand   : {config.reference_ligand}")
    print(f"Candidate ligands  : {len(audit.candidate_ids)}")
    print(f"Validated ligands  : {len(audit.validated_ids)}")
    print(f"Documented exceptions: {len(audit.exception_ids)}")
    print(f"Reported ligands   : {len(audit.report_ids)}")
    print(f"Expected validated : {config.expected_ligands}")
    print(f"Status             : {'PASS' if audit.ok else 'FAIL'}")

    for warning in audit.warnings:
        print(f"WARNING: {warning}")
    for error in audit.errors:
        print(f"ERROR: {error}")

    print("=" * 80)
    return 0 if audit.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
