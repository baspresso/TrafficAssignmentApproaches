"""High-level bilevel Continuous Network Design Problem (CNDP) solving."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._core import (
    BilevelCND,
    LinkConstraint,
    MetricsConfig,
    Network,
    StepConfig,
    TrafficAssignmentApproach,
)
from .datasets import default_constraints_path, load_constraints, load_network
from .tap import make_approach

# Links whose BPR delay cannot respond to capacity changes are pinned to their
# lower bound, mirroring cndp_solver's design-variable filtering.
_INSENSITIVE_EPS = 1e-10


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
               approach: str | TrafficAssignmentApproach = "tapas",
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
    File outputs (trace CSV, metadata JSON, summary CSV) are written only when a
    ``metrics`` config/dict is given; ``final_diagnostics=True`` additionally
    runs the post-run optimality-condition diagnostic and solution CSV dumps.
    The upper-level objective is ``total_travel_time + theta * sum(cap - lb)``
    subject to ``theta * sum(cap - lb) <= budget``.
    """
    steps = [_as_step(entry) for entry in pipeline]
    if not steps:
        raise ValueError("pipeline must contain at least one step")

    if isinstance(network, str):
        network = load_network(network, data_root=data_root)

    if constraints is None:
        constraints = default_constraints_path(network.name, data_root=data_root)
    if isinstance(constraints, (str, Path)):
        constraints_path = Path(constraints)
        if not constraints_path.is_file():
            raise FileNotFoundError(f"Constraints file does not exist: {constraints_path}")
        constraints = load_constraints(constraints_path)
    # Work on copies: the insensitive-link filter must not mutate caller objects.
    constraints = [
        LinkConstraint(c.init_node, c.term_node, c.lower_bound, c.upper_bound,
                       c.investment_cost_param)
        for c in constraints
    ]
    if len(constraints) != network.number_of_links:
        raise ValueError(
            f"constraints ({len(constraints)}) must have one entry per link "
            f"({network.number_of_links})"
        )

    if exclude_insensitive_links:
        for i, constraint in enumerate(constraints):
            link = network.link(i)
            if (abs(link.b) < _INSENSITIVE_EPS or abs(link.power) < _INSENSITIVE_EPS
                    or abs(link.free_flow_time) < _INSENSITIVE_EPS):
                constraint.upper_bound = constraint.lower_bound

    if isinstance(approach, TrafficAssignmentApproach):
        if approach_options:
            raise TypeError("approach_options are only valid when approach is given by name")
        approach_obj = approach
    else:
        approach_obj = make_approach(network, approach, **(approach_options or {}))

    solver = BilevelCND(
        network, approach_obj, constraints, steps,
        link_capacity_selection_threshold=link_threshold,
        budget_threshold=budget_threshold,
        budget_function_multiplier=theta,
        budget_upper_bound=budget,
        metrics=_as_metrics(metrics) if metrics is not None else MetricsConfig(),
        route_search_threads=route_search_threads,
        progress_format=progress,
    )
    solver.set_verbose(verbose)
    solver.set_statistics_enabled(metrics is not None)
    solver.set_final_diagnostics_enabled(final_diagnostics)

    start = time.perf_counter()
    solver.compute_network_design()
    elapsed_seconds = time.perf_counter() - start

    capacities = network.capacities()
    lower_bounds = np.array([c.lower_bound for c in constraints])
    upper_bounds = np.array([c.upper_bound for c in constraints])
    budget_value = theta * float(np.sum(capacities - lower_bounds))
    total_travel_time = network.total_travel_time()
    return CndpResult(
        objective=total_travel_time + budget_value,
        total_travel_time=total_travel_time,
        budget=budget_value,
        budget_upper_bound=budget,
        capacities=capacities,
        flows=network.flows(),
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        elapsed_seconds=elapsed_seconds,
        network=network,
    )
