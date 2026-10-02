"""CSV sample table handling for the Celatlas Python runner."""

from __future__ import annotations

import csv
from pathlib import Path


class SampleTableError(ValueError):
    """Raised when a sample CSV cannot be parsed or validated."""


REQUIRED_COLUMNS = (
    "workflow",
    "chip_number",
    "casno",
    "chemistry",
    "species",
    "method",
    "mode",
)


def read_sample_table(path: str | Path) -> list[dict[str, str]]:
    csv_path = Path(path).expanduser()
    if not csv_path.exists():
        raise SampleTableError(f"Sample CSV not found: {csv_path}")

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise SampleTableError(f"Sample CSV has no header: {csv_path}")
        fieldnames = [_clean_key(name) for name in reader.fieldnames]
        missing = [column for column in REQUIRED_COLUMNS if column not in fieldnames]
        if missing:
            raise SampleTableError(
                f"Sample CSV missing required column(s): {', '.join(missing)}"
            )

        rows: list[dict[str, str]] = []
        for index, raw_row in enumerate(reader, start=2):
            row: dict[str, str] = {"_line_number": str(index)}
            for raw_key, raw_value in raw_row.items():
                if raw_key is None:
                    continue
                key = _clean_key(raw_key)
                row[key] = (raw_value or "").strip()
            if any(value for key, value in row.items() if not key.startswith("_")):
                rows.append(row)
    return rows


def is_enabled(row: dict[str, str]) -> bool:
    value = row.get("enabled", "1").strip().lower()
    return value in {"", "1", "true", "yes", "y", "on", "enable", "enabled"}


def _clean_key(key: str) -> str:
    return key.strip().lstrip("\ufeff").replace("-", "_").lower()
