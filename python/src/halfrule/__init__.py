"""Optimal eligibility criteria for randomised trials."""

from .core import (
    Eligibility,
    design_power,
    eligibility,
    inclusion_bound,
    power_elasticity,
    relative_duration,
    rule_threshold,
    subgroup_threshold,
    time_window,
)

__all__ = [
    "Eligibility",
    "design_power",
    "eligibility",
    "inclusion_bound",
    "power_elasticity",
    "relative_duration",
    "rule_threshold",
    "subgroup_threshold",
    "time_window",
]
__version__ = "0.1.0"
