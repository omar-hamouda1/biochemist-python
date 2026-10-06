"""Validated configuration model for the standardized docking pipeline."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Any, Mapping

import yaml


class DockingConfigError(ValueError):
    """Raised when the docking configuration is missing or invalid."""


@dataclass(frozen=True)
class DockingConfig:
    """Validated standardized docking configuration.

    All project-relative paths are resolved against the repository root, which
    is inferred from the location of ``configs/docking_config.yml``.
    """

    target_pdb: str
    reference_ligand: str
    receptor: Path
    ligand_dir: Path
    output_dir: Path
    manifest: Path
    report: Path
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    exhaustiveness: int
    num_modes: int
    seed: int
    timeout_seconds: int
    smina_executable: str
    expected_ligands: int

    @property
    def center(self) -> tuple[float, float, float]:
        return (self.center_x, self.center_y, self.center_z)

    @property
    def box_size(self) -> tuple[float, float, float]:
        return (self.size_x, self.size_y, self.size_z)


def _require(mapping: Mapping[str, Any], key: str) -> Any:
    value = mapping.get(key)
    if value is None or value == "":
        raise DockingConfigError(f"Missing required docking config key: {key}")
    return value


def _float(mapping: Mapping[str, Any], key: str) -> float:
    try:
        value = float(_require(mapping, key))
    except (TypeError, ValueError) as exc:
        raise DockingConfigError(f"{key} must be a number") from exc
    if not isfinite(value):
        raise DockingConfigError(f"{key} must be finite, got {value}")
    return value


def _positive_float(mapping: Mapping[str, Any], key: str) -> float:
    value = _float(mapping, key)
    if value <= 0:
        raise DockingConfigError(f"{key} must be > 0, got {value}")
    return value


def _positive_int(mapping: Mapping[str, Any], key: str) -> int:
    raw = _require(mapping, key)
    if isinstance(raw, bool):
        raise DockingConfigError(f"{key} must be an integer")
    try:
        value = int(raw)
    except (TypeError, ValueError) as exc:
        raise DockingConfigError(f"{key} must be an integer") from exc
    if isinstance(raw, float) and raw != value:
        raise DockingConfigError(f"{key} must be an integer")
    if value <= 0:
        raise DockingConfigError(f"{key} must be > 0, got {value}")
    return value


def _resolve_project_path(project_root: Path, raw_path: str) -> Path:
    path = Path(raw_path).expanduser()
    return path if path.is_absolute() else project_root / path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "configs" / "docking_config.yml"


def load_docking_config(config_path: str | Path = DEFAULT_CONFIG_PATH) -> DockingConfig:
    """Load and validate the authoritative standardized docking configuration."""

    path = Path(config_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Docking config not found: {path}")

    with path.open() as handle:
        raw = yaml.safe_load(handle) or {}

    if not isinstance(raw, dict):
        raise DockingConfigError("Docking config must contain a YAML mapping")

    project_root = path.parent.parent
    executable = str(raw.get("smina_executable", "smina")).strip()
    if not executable:
        raise DockingConfigError("smina_executable must not be empty")

    expected_ligands = _positive_int(raw, "expected_ligands")

    return DockingConfig(
        target_pdb=str(_require(raw, "target_pdb")),
        reference_ligand=str(_require(raw, "reference_ligand")),
        receptor=_resolve_project_path(project_root, str(_require(raw, "receptor"))),
        ligand_dir=_resolve_project_path(project_root, str(_require(raw, "ligand_dir"))),
        output_dir=_resolve_project_path(project_root, str(_require(raw, "output_dir"))),
        manifest=_resolve_project_path(project_root, str(_require(raw, "manifest"))),
        report=_resolve_project_path(project_root, str(_require(raw, "report"))),
        center_x=_float(raw, "center_x"),
        center_y=_float(raw, "center_y"),
        center_z=_float(raw, "center_z"),
        size_x=_positive_float(raw, "size_x"),
        size_y=_positive_float(raw, "size_y"),
        size_z=_positive_float(raw, "size_z"),
        exhaustiveness=_positive_int(raw, "exhaustiveness"),
        num_modes=_positive_int(raw, "num_modes"),
        seed=_positive_int(raw, "seed"),
        timeout_seconds=_positive_int(raw, "timeout_seconds"),
        smina_executable=executable,
        expected_ligands=expected_ligands,
    )
