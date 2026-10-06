"""Provenance helpers for reproducible structure-preparation and docking steps."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def executable_version(executable: str) -> str:
    """Return normalized executable version text."""
    result = subprocess.run(
        [executable, "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    return " ".join((result.stdout or result.stderr).strip().split())


def provenance_path(output_path: Path) -> Path:
    """Return the JSON sidecar path associated with an output artifact."""
    return output_path.with_name(f"{output_path.name}.provenance.json")


def _fingerprint_payload(metadata: dict[str, Any]) -> dict[str, Any]:
    """Return fields that define the computation independently of output bytes."""
    return {
        "schema_version": metadata["schema_version"],
        "source_sha256": metadata["source_sha256"],
        "producer_sha256": metadata["producer_sha256"],
        "tool_versions": metadata["tool_versions"],
        "protocol": metadata["protocol"],
    }


def build_fingerprint(metadata: dict[str, Any]) -> str:
    """Build a stable computation fingerprint."""
    encoded = json.dumps(
        _fingerprint_payload(metadata),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_artifact_metadata(
    *,
    source_path: Path,
    producer_path: Path,
    tool_versions: dict[str, str],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    """Build metadata for a derived structural artifact."""
    metadata: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
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
    output_path: Path,
    affinity: float | None = None,
) -> None:
    """Write metadata including the exact output-file hash."""
    payload = {
        **metadata,
        "output_sha256": sha256_file(output_path),
    }
    if affinity is not None:
        payload["affinity"] = affinity
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "
",
        encoding="utf-8",
    )


def validate_existing_artifact(
    *,
    output_path: Path,
    metadata_path: Path,
    expected_metadata: dict[str, Any],
) -> bool:
    """Return True only when output bytes and computation fingerprint both match."""
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
    """Backward-compatible docking-output validation including affinity."""
    if not validate_existing_artifact(
        output_path=output_path,
        metadata_path=metadata_path,
        expected_metadata=expected_metadata,
    ):
        return False

    metadata = load_metadata(metadata_path)
    if expected_affinity is None:
        return True
    try:
        return float(metadata.get("affinity")) == expected_affinity
    except (TypeError, ValueError):
        return False


def output_has_atoms(path: Path) -> bool:
    """Return True when a PDBQT contains at least one atom record."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        return any(line.startswith(("ATOM", "HETATM")) for line in handle)
