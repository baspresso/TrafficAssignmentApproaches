"""High-level bilevel Continuous Network Design Problem (CNDP) solving."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._core import (
    BilevelCND,
    CndpOptions,
    MetricsConfig,
    Network,
    StepConfig,
    TapOptions,
    TrafficAssignmentApproach,
)
from .datasets import _resolve_constraints, _resolve_network
from .tap import make_approach


def step(type: str, **options) -> StepConfig:
    """Build a pipeline step, e.g. ``step("nlopt", algorithm="LN_COBYLA")``.

    Option names match the ``[[pipeline]]`` TOML keys (``max_iterations``,
    ``tolerance``, ``algorithm``, ``population_size``, ...).
    """
    config = StepConfig()
    config.type = type
    for key, value in options.items():
        if not hasattr(config, key):
            raise ValueError(f"Unknown pipeline step option '{key}'")
        setattr(config, key, value)
    return config


def _as_step(value) -> StepConfig:
    if isinstance(value, StepConfig):
        return value
    if isinstance(value, dict):
        options = dict(value)
        try:
            step_type = options.pop("type")
        except KeyError:
            raise ValueError("pipeline step dicts need a 'type' key") from None
        return step(step_type, **options)
    raise TypeError(f"pipeline entries must be StepConfig or dict, got {value!r}")


def _as_metrics(metrics) -> MetricsConfig:
    if isinstance(metrics, MetricsConfig):
        return metrics
    config = MetricsConfig()
    for key, value in dict(metrics).items():
        if not hasattr(config, key):
            raise ValueError(f"Unknown metrics option '{key}'")
        setattr(config, key, value)
    return config


@dataclass
class CndpResult:
    """Outcome of a bilevel CNDP solve (final user-equilibrium state)."""

    objective: float
    total_travel_time: float
    budget: float
    budget_upper_bound: float
    capacities: np.ndarray
    flows: np.ndarray
    lower_bounds: np.ndarray
    upper_bounds: np.ndarray
    elapsed_seconds: float
    network: Network


def solve_cndp(network: Network | str, pipeline, constraints=None, *,
               approach: str | TapOptions | TrafficAssignmentApproach = "tapas",
               approach_options: dict | None = None, theta: float = 5.0,
               budget: float = 100000.0, budget_threshold: float = 0.1,
               link_threshold: float = 1e-3, route_search_threads: int = 1,
               metrics=None, progress: str = "none", verbose: bool = False,
               exclude_insensitive_links: bool = True,
               final_diagnostics: bool = False, data_root=None) -> CndpResult:
    """Solve the bilevel CNDP; defaults mirror ``cndp_solver`` with file output off.

    ``pipeline`` is a non-empty sequence of :func:`step` results or dicts with a
    ``"type"`` key. ``constraints`` is a list of :class:`LinkConstraint`, a CSV
    path, or ``None`` to use the dataset's conventional constraints file.
    File outputs (trace CSV, metadata JSON, summary and final route-count CSVs)
    are written only when a ``metrics`` config/dict is given;
    ``final_diagnostics=True`` additionally
    runs the post-run optimality-condition diagnostic and solution CSV dumps.
    The upper-level objective is ``total_travel_time + theta * sum(cap - lb)``
    subject to ``theta * sum(cap - lb) <= budget``.
    """
    steps = [_as_step(entry) for entry in pipeline]
    if not steps:
        raise ValueError("pipeline must contain at least one step")

    network = _resolve_network(network, data_root=data_root)
    constraints = _resolve_constraints(network, constraints, data_root=data_root)
    if isinstance(approach, TrafficAssignmentApproach):
        if approach_options:
            raise TypeError("approach_options are only valid when approach is given by name")
        approach_obj = approach
    else:
        approach_obj = make_approach(network, approach, **(approach_options or {}))

    options = CndpOptions()
    options.link_capacity_selection_threshold = link_threshold
    options.budget_threshold = budget_threshold
    options.budget_function_multiplier = theta
    options.budget_upper_bound = budget
    options.route_search_threads = route_search_threads
    options.progress_format = progress
    options.exclude_insensitive_links = exclude_insensitive_links
    options.statistics_enabled = metrics is not None
    options.final_diagnostics = final_diagnostics
    options.verbose = verbose
    if metrics is not None:
        options.metrics = _as_metrics(metrics)

    solver = BilevelCND(network, approach_obj, constraints, steps, options)
    native = solver.compute_network_design()
    return CndpResult(
        objective=native.objective,
        total_travel_time=native.total_travel_time,
        budget=native.budget,
        budget_upper_bound=native.budget_upper_bound,
        capacities=native.capacities,
        flows=native.flows,
        lower_bounds=native.lower_bounds,
        upper_bounds=native.upper_bounds,
        elapsed_seconds=native.elapsed_seconds,
        network=network,
    )
