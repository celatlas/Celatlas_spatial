"""Configuration loading helpers for the Celatlas Python runner."""

from __future__ import annotations

import json
import os
import re
import shlex
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when a runner config cannot be loaded."""


ENV_TO_CONFIG_KEY = {
    "CELATLAS_WORKSPACE": "workspace",
    "CELATLAS_RESULTS_ROOT": "results_root",
    "CELATLAS_FASTQ_ROOT": "fastq_root",
    "CELATLAS_MASK_DIR": "mask_dir",
    "CELATLAS_IMAGE_DIR": "image_dir",
    "CELATLAS_REFERENCE_DIR": "reference_dir",
    "CELATLAS_SX_REFERENCE_DIR": "sx_reference_dir",
    "CELATLAS_FFPE_REFERENCE_DIR": "ffpe_reference_dir",
    "CELATLAS_SRC_DIR": "src_dir",
    "CELATLAS_ENV_NAME": "env_name",
    "CELATLAS_DEFAULT_THREAD": "default_thread",
    "CELATLAS_DEFAULT_BIN": "default_bin",
    "CELATLAS_DEFAULT_PIXEL_SIZE": "default_pixel_size",
}

CONFIG_TO_ENV = {value: key for key, value in ENV_TO_CONFIG_KEY.items()}


def load_config(path: str | os.PathLike[str] | None) -> dict[str, Any]:
    """Load a YAML, JSON, or simple shell env config file."""

    if not path:
        return {}

    config_path = Path(path).expanduser()
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")

    suffix = config_path.suffix.lower()
    name = config_path.name.lower()
    if suffix in {".yaml", ".yml"} or name.endswith((".yaml.example", ".yml.example")):
        return _load_yaml(config_path)
    if suffix == ".json" or name.endswith(".json.example"):
        return _load_json(config_path)
    return _load_env_config(config_path)


def load_profile(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Load a profile file into config/env values."""

    return load_config(path)


def config_to_env(config: dict[str, Any]) -> dict[str, str]:
    """Convert runner config keys and explicit env values into environment variables."""

    env: dict[str, str] = {}

    for key, value in (config.get("env") or {}).items():
        if value is not None and str(value) != "":
            env[str(key)] = str(value)

    for key, env_key in CONFIG_TO_ENV.items():
        value = get_config_value(config, key)
        if value is not None and str(value) != "":
            env[env_key] = str(value)

    for key, value in config.items():
        if key.startswith("CELATLAS_") and value is not None and str(value) != "":
            env[key] = str(value)

    return env


def get_config_value(config: dict[str, Any], *keys: str) -> Any:
    """Return the first present value from flat config or common nested sections."""

    sections = (
        "defaults",
        "paths",
        "runner",
        "env",
    )
    for key in keys:
        if key in config and config[key] not in (None, ""):
            return config[key]
        env_key = CONFIG_TO_ENV.get(key)
        if env_key and env_key in config and config[env_key] not in (None, ""):
            return config[env_key]
        for section in sections:
            data = config.get(section)
            if isinstance(data, dict):
                if key in data and data[key] not in (None, ""):
                    return data[key]
                if env_key and env_key in data and data[env_key] not in (None, ""):
                    return data[env_key]
    return None


def merge_configs(*configs: dict[str, Any]) -> dict[str, Any]:
    """Merge config dictionaries. Later values override earlier values."""

    merged: dict[str, Any] = {}
    for config in configs:
        for key, value in config.items():
            if isinstance(value, dict) and isinstance(merged.get(key), dict):
                nested = dict(merged[key])
                nested.update(value)
                merged[key] = nested
            elif value is not None and value != "":
                merged[key] = value
    return merged


def _load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise ConfigError("PyYAML is required to read YAML runner configs") from exc

    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"YAML config must contain a mapping: {path}")
    return data


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ConfigError(f"JSON config must contain an object: {path}")
    return data


def _load_env_config(path: Path) -> dict[str, Any]:
    env_values = read_env_file(path)
    config: dict[str, Any] = {"env": env_values}
    for env_key, config_key in ENV_TO_CONFIG_KEY.items():
        if env_key in env_values:
            config[config_key] = env_values[env_key]
    return config


def read_env_file(path: str | os.PathLike[str]) -> dict[str, str]:
    """Read simple KEY=value or export KEY=value lines.

    This parser intentionally supports the subset used by Celatlas env/profile
    files. It ignores shell control-flow lines and does not execute code.
    """

    env_path = Path(path).expanduser()
    values: dict[str, str] = {}
    pattern = re.compile(r"^(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=(.*)$")

    with env_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            match = pattern.match(stripped)
            if not match:
                continue
            key, raw_value = match.groups()
            value = _parse_shell_value(raw_value)
            values[key] = _expand_shell_vars(value, values)
    return values


def _parse_shell_value(raw_value: str) -> str:
    try:
        parts = shlex.split(raw_value, comments=True, posix=True)
    except ValueError:
        return raw_value.strip().strip("'\"")
    if not parts:
        return ""
    return " ".join(parts)


def _expand_shell_vars(value: str, local_values: dict[str, str]) -> str:
    def replace_default(match: re.Match[str]) -> str:
        key = match.group(1)
        default = match.group(2)
        current = local_values.get(key, os.environ.get(key, ""))
        if current:
            return current
        return _expand_shell_vars(default or "", local_values)

    default_pattern = re.compile(
        r"\$\{([A-Za-z_][A-Za-z0-9_]*):-((?:[^{}]|\$\{[^{}]+\})*)\}"
    )
    previous = None
    while previous != value:
        previous = value
        value = default_pattern.sub(replace_default, value)

    def replace_plain(match: re.Match[str]) -> str:
        key = match.group(1)
        return local_values.get(key, os.environ.get(key, ""))

    value = re.sub(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}", replace_plain, value)
    value = re.sub(r"\$([A-Za-z_][A-Za-z0-9_]*)", replace_plain, value)
    return value
