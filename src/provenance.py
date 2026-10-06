"""Provenance helpers for reproducible structure preparation and docking."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 2


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def executable_version(executable: str) -> str:
    """Return normalized executable version text, or unavailable."""
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except (FileNotFoundError, OSError):
        return "unavailable"

    output = (result.stdout or result.stderr or "").strip()
    return " ".join(output.split()) or "unavailable"


def provenance_path(output_path: Path) -> Path:
    """Return the JSON sidecar path associated with an output artifact."""
    return output_path.with_name(f"{output_path.name}.provenance.json")


def _fingerprint_payload(metadata: dict[str, Any]) -> dict[str, Any]:
    """Return path-independent fields that define the computation."""
    keys = (
        "schema_version",
        "artifact_kind",
        "config_sha256",
        "source_sha256",
        "receptor_sha256",
        "ligand_sha256",
        "producer_sha256",
        "runner_sha256",
        "parser_sha256",
        "smina_version",
        "tool_versions",
        "protocol",
    )

    payload = {key: metadata[key] for key in keys if key in metadata}

    protocol = payload.get("protocol")
    if isinstance(protocol, dict):
        protocol = dict(protocol)
        protocol.pop("command", None)
        payload["protocol"] = protocol

    return payload


def build_fingerprint(metadata: dict[str, Any]) -> str:
    """Build a stable path-independent computation fingerprint."""
    encoded = json.dumps(
        _fingerprint_payload(metadata),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    return hashlib.sha256(encoded).hexdigest()


def build_metadata(
    *,
    config_path: Path,
    receptor_path: Path,
    ligand_path: Path,
    runner_path: Path,
    parser_path: Path,
    smina_version: str,
    command: list[str],
) -> dict[str, Any]:
    """Build provenance metadata for one standardized docking pose."""

    metadata: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "artifact_kind": "standardized_docking_pose",
        "config_sha256": sha256_file(config_path),
        "receptor_sha256": sha256_file(receptor_path),
        "ligand_sha256": sha256_file(ligand_path),
        "runner_sha256": sha256_file(runner_path),
        "parser_sha256": sha256_file(parser_path),
        "smina_version": smina_version,
        "protocol": {
            "command": list(command),
            "command_paths_are_excluded_from_fingerprint": True,
        },
    }

    metadata["fingerprint"] = build_fingerprint(metadata)

    return metadata


def build_artifact_metadata(
    *,
    source_path: Path,
    producer_path: Path,
    tool_versions: dict[str, str],
    protocol: dict[str, Any],
    artifact_kind: str,
) -> dict[str, Any]:
    """Build metadata for a generic derived structural artifact."""

    metadata: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "artifact_kind": artifact_kind,
        "source_sha256": sha256_file(source_path),
        "producer_sha256": sha256_file(producer_path),
        "tool_versions": dict(sorted(tool_versions.items())),
        "protocol": protocol,
    }

    metadata["fingerprint"] = build_fingerprint(metadata)

    return metadata


def load_metadata(path: Path) -> dict[str, Any] | None:
    """Load a valid JSON sidecar, returning None when missing or malformed."""

    if not path.exists() or path.stat().st_size == 0:
        return None

    try:
        with path.open(encoding="utf-8") as handle:
            metadata = json.load(handle)

    except (OSError, json.JSONDecodeError):
        return None

    return metadata if isinstance(metadata, dict) else None


def write_metadata(
    path: Path,
    metadata: dict[str, Any],
    *,
    output_path: Path | None = None,
    output_sha256: str | None = None,
    affinity: float | None = None,
    extra: dict[str, Any] | None = None,
) -> None:
    """Write metadata including the exact output-file hash."""

    if output_sha256 is None:
        if output_path is None:
            raise ValueError("output_path or output_sha256 is required")

        output_sha256 = sha256_file(output_path)

    payload: dict[str, Any] = {
        **metadata,
        "output_sha256": output_sha256,
    }

    if affinity is not None:
        payload["affinity"] = affinity

    if extra:
        payload.update(extra)

    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def validate_existing_artifact(
    *,
    output_path: Path,
    metadata_path: Path,
    expected_metadata: dict[str, Any],
) -> bool:
    """Return True only when fingerprint and output bytes both match."""

    metadata = load_metadata(metadata_path)

    if metadata is None:
        return False

    if metadata.get("fingerprint") != expected_metadata.get("fingerprint"):
        return False

    return metadata.get("output_sha256") == sha256_file(output_path)


def validate_existing_output(
    *,
    output_path: Path,
    metadata_path: Path,
    expected_metadata: dict[str, Any],
    expected_affinity: float | None = None,
) -> bool:
    """Validate an existing docking output, including its recorded affinity."""

    if not validate_existing_artifact(
        output_path=output_path,
        metadata_path=metadata_path,
        expected_metadata=expected_metadata,
    ):
        return False

    if expected_affinity is None:
        return True

    metadata = load_metadata(metadata_path)

    if metadata is None:
        return False

    try:
        recorded_affinity = metadata.get("affinity")

        if recorded_affinity is None:
            return False

        return float(recorded_affinity) == expected_affinity

    except (TypeError, ValueError):
        return False


def output_has_atoms(path: Path) -> bool:
    """Return True when a PDBQT contains at least one atom record."""

    with path.open(encoding="utf-8", errors="replace") as handle:
        return any(
            line.startswith(("ATOM", "HETATM"))
            for line in handle
        )