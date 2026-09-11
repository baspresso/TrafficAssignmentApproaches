"""Solve an analytic example without downloading a dataset.

Run after installing the package: python examples/two_routes.py
"""

import traffic_assignment as ta

network = ta.network_from_arrays(
    "TwoRoutes", init_node=[0, 0, 2], term_node=[1, 2, 1],
    capacity=10, free_flow_time=[2, 1, 1], b=1, power=1,
    demand=[[0, 10], [0, 0]], n_nodes=3,
)
options = ta.TapOptions()
options.relative_gap_tolerance = 1e-12
result = ta.solve_tap(network, approach=options)

print("Link flows:", result.flows)  # [5, 5, 5]
print("Total travel time:", result.total_travel_time)  # 30
print("Relative gap:", result.relative_gap)
