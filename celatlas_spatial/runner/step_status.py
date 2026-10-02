"""Step output status helpers for Python runner dry-run/resume planning."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .steps import StepPlan


@dataclass(frozen=True)
class StepOutputStatus:
    path: Path
    exists: bool


@dataclass(frozen=True)
class StepStatus:
    name: str
    state: str
    outputs: tuple[StepOutputStatus, ...]


def evaluate_step_status(plans: list[StepPlan]) -> list[StepStatus]:
    return [_evaluate_one(plan) for plan in plans]


def _evaluate_one(plan: StepPlan) -> StepStatus:
    outputs = tuple(StepOutputStatus(path, path.exists()) for path in plan.outputs)
    if not outputs:
        state = "unknown"
    elif all(output.exists for output in outputs):
        state = "complete"
    elif any(output.exists for output in outputs):
        state = "partial"
    else:
        state = "missing"
    return StepStatus(plan.name, state, outputs)
