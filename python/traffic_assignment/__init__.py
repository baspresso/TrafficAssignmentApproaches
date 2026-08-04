"""Traffic Assignment (TAP) and Continuous Network Design (CNDP) solvers.

Python bindings for the C++ solvers in this repository: route-based and
bush-based (TAPAS) user-equilibrium assignment, and the bilevel CNDP solver
with NLopt / optimality-condition / OptimLib pipeline steps.

Quick start::

    import traffic_assignment as ta

    result = ta.solve_tap("SiouxFalls")
    design = ta.solve_cndp(
        "SiouxFalls",
        [ta.step("nlopt", algorithm="LN_COBYLA", max_iterations=100)],
        budget=10000.0,
    )
"""

from ._core import (
    BilevelCND,
    Link,
    LinkConstraint,
    MetricsConfig,
    Network,
    RouteBasedApproach,
    StepConfig,
    TapasApproach,
    TrafficAssignmentApproach,
)
from .cndp import CndpResult, solve_cndp, step
from .datasets import (
    constraints_from_arrays,
    default_constraints_path,
    find_data_root,
    load_constraints,
    load_network,
    network_from_arrays,
)
from .tap import TapResult, make_approach, solve_tap

__version__ = "0.1.0"

__all__ = [
    "BilevelCND",
    "CndpResult",
    "Link",
    "LinkConstraint",
    "MetricsConfig",
    "Network",
    "RouteBasedApproach",
    "StepConfig",
    "TapResult",
    "TapasApproach",
    "TrafficAssignmentApproach",
    "constraints_from_arrays",
    "default_constraints_path",
    "find_data_root",
    "load_constraints",
    "load_network",
    "make_approach",
    "network_from_arrays",
    "solve_cndp",
    "solve_tap",
    "step",
    "__version__",
]
