"""Capture the exact local Conda and pip environment used for validation."""

from __future__ import annotations

import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT_ROOT / "docking" / "results" / "standardized_environment.txt"


def run(command: list[str]) -> str:
    """Run a command and return its combined output."""
    try:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        return f"COMMAND UNAVAILABLE: {exc}\n"

    body = result.stdout or ""
    if result.stderr:
        body += "\n[stderr]\n" + result.stderr
    return body.strip() + f"\n[exit_code={result.returncode}]\n"


def main() -> None:
    """Write an explicit environment snapshot for reproducibility."""
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(
        [
            "# Exact validation environment snapshot",
            "# Generated locally; do not edit manually.",
            "",
            "===== conda list --explicit =====",
            run(["conda", "list", "--explicit"]),
            "===== python --version =====",
            run(["python", "--version"]),
            "===== pip freeze =====",
            run(["python", "-m", "pip", "freeze"]),
            "===== external tools =====",
            "----- smina -----",
            run(["smina", "--version"]),
            "----- obabel -----",
            run(["obabel", "--version"]),
            "----- pdb2pqr -----",
            run(["pdb2pqr", "--version"]),
        ]
    )
    OUTPUT.write_text(text + "\n", encoding="utf-8")
    print(f"Wrote environment snapshot: {OUTPUT}")


if __name__ == "__main__":
    main()
