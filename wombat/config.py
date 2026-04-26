"""Configuration loader for ConvergentDecidua.

Usage:
    from wombat.config import load_config
    datasets = load_config("datasets")
"""

from __future__ import annotations

from pathlib import Path

import yaml

_CONFIGS_DIR = Path(__file__).resolve().parent.parent / "configs"

_REQUIRED_KEYS: dict[str, list[str]] = {
    "datasets": ["accession", "species", "assay"],
    "species": ["name", "taxon_id"],
    "markers": ["cell_type_markers"],
}


def load_config(name: str) -> dict:
    """Load and validate a YAML config by name.

    Parameters
    ----------
    name : str
        Config name without extension (e.g. ``"datasets"``).

    Returns
    -------
    dict
        Parsed YAML content.

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    ValueError
        If required keys are missing.
    """
    path = _CONFIGS_DIR / f"{name}.yaml"
    if not path.exists():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    with open(path) as fh:
        data = yaml.safe_load(fh)

    if data is None:
        msg = f"Config file is empty: {path}"
        raise ValueError(msg)

    _validate(name, data)
    return data


def _validate(name: str, data: dict | list) -> None:
    """Check that required keys are present."""
    required = _REQUIRED_KEYS.get(name)
    if required is None:
        return

    if isinstance(data, list):
        for i, entry in enumerate(data):
            missing = [k for k in required if k not in entry]
            if missing:
                msg = f"{name}[{i}] missing required keys: {missing}"
                raise ValueError(msg)
    elif isinstance(data, dict):
        # For dict configs, check top-level keys
        missing = [k for k in required if k not in data]
        if missing:
            msg = f"{name} missing required keys: {missing}"
            raise ValueError(msg)


def validate_all() -> list[str]:
    """Validate all config files found in the configs directory.

    Returns
    -------
    list[str]
        List of error messages (empty if all valid).
    """
    errors: list[str] = []
    for path in sorted(_CONFIGS_DIR.glob("*.yaml")):
        name = path.stem
        try:
            load_config(name)
        except (FileNotFoundError, ValueError, yaml.YAMLError) as exc:
            errors.append(f"{name}: {exc}")
    return errors
