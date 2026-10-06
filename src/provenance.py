"""Reproducibility metadata for the standardized docking workflow."""

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
    text = (result.stdout or result.stderr).strip()
    return " ".join(text.split())


def provenance_path(output_path: Path) -> Path:
    """Return the sidecar path associated with a docking output."""
    return output_path.with_name(f"{output_path.name}.provenance.json")


def _fingerprint_payload(metadata: dict[str, Any]) -> dict[str, Any]:
    return {
        key: metadata[key]
        for key in (
            "schema_version",
            "config_sha256",
            "receptor_sha256",
            "ligand_pdbqt_sha256",
            "runner_sha256",
            "parser_sha256",
            "smina_version",
        )
    }


def build_fingerprint(metadata: dict[str, Any]) -> str:
    """Build a stable fingerprint from all computation-defining inputs."""
    payload = _fingerprint_payload(metadata)
    encoded = json.dumps(
        payload,
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
    """Build computation-defining metadata for one docking run."""
    metadata: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "config_sha256": sha256_file(config_path),
        "receptor_sha256": sha256_file(receptor_path),
        "ligand_pdbqt_sha256": sha256_file(ligand_path),
        "runner_sha256": sha256_file(runner_path),
        "parser_sha256": sha256_file(parser_path),
        "smina_version": smina_version,
        "command": command,
    }
    metadata["fingerprint"] = build_fingerprint(metadata)
    return metadata


def write_metadata(path: Path, metadata: dict[str, Any], *, output_sha256: str, affinity: float) -> None:
    """Write validated output metadata as a JSON sidecar."""
    payload = {
        **metadata,
        "output_sha256": output_sha256,
        "affinity": affinity,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_metadata(path: Path) -> dict[str, Any] | None:
    """Load a provenance sidecar, returning None when invalid or missing."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        with path.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return metadata if isinstance(metadata, dict) else None


def validate_existing_output(
    *,
    output_path: Path,
    metadata_path: Path,
    expected_metadata: dict[str, Any],
    expected_affinity: float | None = None,
) -> bool:
    """Return True only when output and provenance match current inputs."""
    metadata = load_metadata(metadata_path)
    if metadata is None:
        return False

    if metadata.get("fingerprint") != expected_metadata.get("fingerprint"):
        return False

    if metadata.get("output_sha256") != sha256_file(output_path):
        return False

    if expected_affinity is not None:
        try:
            if float(metadata.get("affinity")) != expected_affinity:
                return False
        except (TypeError, ValueError):
            return False

    return True


def output_has_atoms(path: Path) -> bool:
    """Return True when a PDBQT contains at least one atom record."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        return any(
            line.startswith(("ATOM", "HETATM"))
            for line in handle
        )
