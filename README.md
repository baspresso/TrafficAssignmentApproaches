# TrafficAssignmentApproaches

This repository presents a set of algorithmic implementations designed to solve the classical **Traffic Assignment Problem (TAP)**. The work stems from graduate-level research conducted at Saint Petersburg State University and includes both standard and novel approaches to equilibrium-based flow distribution on transportation networks.

## Problem Overview

The Traffic Assignment Problem (TAP) involves determining an equilibrium distribution of traffic flows across a transportation network, such that no traveler can reduce their travel time by unilaterally changing routes. This condition corresponds to **Wardrop's first principle of user equilibrium**. TAP is fundamental to transportation systems analysis and infrastructure planning, as it captures the interaction between route choices and network congestion.

Mathematically, TAP is often formulated as a convex optimization problem with nonlinear cost functions, representing travel time as a function of link load. Solving this problem with high precision is particularly important in bilevel applications such as network design and congestion pricing.

## Implemented Approaches

The repository includes the following classes of traffic assignment methods:

- **Path-based Methods**
  Manage flow distributions across explicit paths for each origin–destination (OD) pair. This approach supports high-precision computation and direct access to route-level flow information. The implementation follows a route-equilibration strategy grounded in recent academic literature.


- **TAPAS (Traffic Assignment by Paired Alternative Segments)**
  A refinement of bush-based algorithms that employs paired segments with common endpoints but disjoint paths. TAPAS supports both line-search and Newton-step flow updates, allowing for fast and accurate convergence. This method is especially effective in high-precision applications.

- **Demand-based Method (Novel Contribution)**
  A newly developed algorithm introduced by the author, which models dynamic redistribution of OD demand through the network. The approach uses a pump-out mechanism to iteratively reassign unmet demand. While preliminary results indicate lower accuracy and slower convergence, the method presents a promising direction for further investigation.

## Technical Implementation

- **Programming Languages**
  - Core algorithm implementations: **C++20**
  - Data preprocessing and experiment orchestration: **Python**

- **Build System**
  CMake (>= 3.14) with Ninja, using system-installed Eigen3, NLopt, and Boost.

## Build and Run

### Prerequisites (Ubuntu/Debian)

```bash
sudo apt install build-essential cmake ninja-build pkg-config \
                 libeigen3-dev libnlopt-cxx-dev libboost-all-dev
```

(`toml++` and `OptimLib` are fetched automatically by CMake `FetchContent`.)

### Configure

```bash
cmake --preset linux-release           # or linux-debug, linux-relwithdebinfo
```

### Build

```bash
cmake --build --preset build-release   # or build-debug, build-relwithdebinfo
```

### Run

```bash
./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
./build/linux-release/cndp_solver --help
```

A separate standalone TAP solver is also produced:

```bash
./build/linux-release/tap_solver --config configs/tap.siouxfalls.toml
```

## Layered Runtime Config (defaults → config → env → CLI)

`cndp_solver` supports layered configuration:

- Defaults in code
- TOML config file (`--config ...` or `CND_CONFIG`)
- Environment variable overrides (`CND_*`)
- CLI overrides (highest priority)

Example config file: [`configs/cnd.siouxfalls.toml`](configs/cnd.siouxfalls.toml).

Run using config only:

```bash
./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
```

Run with environment override:

```bash
CND_ROUTE_THREADS=1 ./build/linux-release/cndp_solver --config configs/cnd.siouxfalls.toml
```

Run with CLI override (wins over env + file):

```bash
CND_ROUTE_THREADS=1 ./build/linux-release/cndp_solver \
    --config configs/cnd.siouxfalls.toml \
    --route-threads 2 \
    --max-standard-iters 150
```

Show all available options:

```bash
./build/linux-release/cndp_solver --help
```

## Benchmark Networks

The following standardized test networks are used to evaluate algorithmic performance:

- **Sioux Falls Network**
  24 nodes, 24 zones, 76 links

- **Anaheim Network**
  416 nodes, 38 zones, 914 links

These networks are sourced from the open-access [Transportation Networks for Research](https://github.com/bstabler/TransportationNetworks) repository and conform to the TNTP format. Each instance includes a demand matrix and BPR-formulated travel time functions.

## Experimental Results (Summary)

The implemented algorithms were evaluated with respect to solution precision (measured by the Relative Gap, RGAP) and computation time. Key findings include:

- The **TAPAS-NewtonStep** method achieved the best performance, reaching machine-precision accuracy (RGAP ≈ 10⁻¹⁵) in the shortest runtime.
- The **TAPAS-LineSearch** variant also performed favorably, with slightly longer convergence times but similarly high precision.
- The **Path-based** implementation achieved acceptable precision (RGAP ≈ 10⁻¹²), though it incurred higher computational cost due to the maintenance of explicit paths.
- The **Demand-based** method, in its current form, did not achieve comparable precision and exhibited slower convergence, suggesting it may benefit from further optimization and theoretical refinement.
