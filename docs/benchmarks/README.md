# Traffic assignment benchmark study

Measured 2026-10-13 · **full profile** · 78/78 measured runs passed validation.

This study compares the repository's implementations under identical inputs and requested accuracy. Results apply to this machine and build. External solver comparisons and global CNDP optimality are outside this study.

## Inputs and model

| Dataset | Nodes | Links | Interzonal OD pairs | Source demand | Road-network demand |
| --- | ---: | ---: | ---: | ---: | ---: |
| SiouxFalls | 24 | 76 | 528 | 360,600.00 | 360,600.00 |
| ChicagoSketch | 933 | 2950 | 93135 | 1,260,907.44 | 1,137,493.44 |

Source: [Transportation Networks for Research](https://github.com/bstabler/TransportationNetworks/tree/375e0da93858c547230c5cf9ea8a96de4ccff29e). Exact revision and input checksums are in [manifest.json](manifest.json). Fractional demand is preserved. Intrazonal trips remain in the input as zero-link trips; they contribute no road flow or travel time. Road-network demand excludes them.

The selected networks use BPR travel-time costs and allow transit through every node. **Chicago Sketch is a travel-time-only variant**: toll and distance terms from the published generalized-cost reference are excluded, so its published reference objective does not apply here. Sioux Falls is a small research/debugging benchmark. Travel-time totals retain native dataset units.

## Assignment accuracy and runtime

Times are medians with 25th–75th percentiles, in seconds. Statistics use successful measured runs; the count exposes incomplete configurations. Warm-ups are excluded.

| Dataset | Method | Target gap | Successful runs | Median [Q25, Q75] (s) | Max achieved gap |
| --- | --- | ---: | ---: | ---: | ---: |
| SiouxFalls | tapas | 1e-06 | 5/5 | 0.0144 [0.0143, 0.0255] | 6.41e-07 |
| SiouxFalls | routebased | 1e-06 | 5/5 | 0.0131 [0.0122, 0.0187] | 9.07e-07 |
| SiouxFalls | tapas | 1e-08 | 5/5 | 0.0410 [0.0370, 0.0528] | 2.85e-10 |
| SiouxFalls | routebased | 1e-08 | 5/5 | 0.0444 [0.0438, 0.0451] | 4.11e-09 |
| SiouxFalls | tapas | 1e-10 | 5/5 | 0.0594 [0.0571, 0.0617] | 1.67e-12 |
| SiouxFalls | routebased | 1e-10 | 5/5 | 0.0660 [0.0551, 0.0662] | 6.42e-11 |
| ChicagoSketch | tapas | 1e-06 | 5/5 | 3.6315 [3.6203, 3.6516] | 1.70e-07 |
| ChicagoSketch | routebased | 1e-06 | 5/5 | 2.0868 [2.0404, 2.1114] | 9.97e-08 |
| ChicagoSketch | tapas | 1e-08 | 5/5 | 4.8484 [4.8386, 4.8936] | 2.76e-18 |
| ChicagoSketch | routebased | 1e-08 | 5/5 | 3.9594 [3.8113, 4.0087] | 2.75e-10 |
| ChicagoSketch | tapas | 1e-10 | 5/5 | 4.9105 [4.8251, 4.9153] | 2.76e-18 |
| ChicagoSketch | routebased | 1e-10 | 5/5 | 5.8035 [5.6894, 5.8498] | 6.13e-12 |

![Runtime at matched accuracy](figures/assignment_runtime.svg)

## Parallel route-search scaling

Chicago Sketch, route-based solver, target relative gap 1e-8. Speedup measures the entire solve: flow redistribution remains sequential. Speedups below 1 indicate a slowdown. Only fully validated configurations with matching travel-time and Beckmann objectives receive a speedup.

| Threads | Successful runs | Median [Q25, Q75] (s) | Speedup |
| ---: | ---: | ---: | ---: |
| 1 | 5/5 | 3.9594 [3.8113, 4.0087] | 1.00× |
| 2 | 5/5 | 4.4194 [4.3358, 4.4684] | 0.90× |
| 4 | 5/5 | 4.1089 [4.0583, 4.1530] | 0.96× |
| 8 | 5/5 | 3.8965 [3.8678, 4.0174] | 1.02× |

![Parallel scaling](figures/parallel_scaling.svg)

## Capacity-optimization case

Sioux Falls with COBYLA, at most 500 objective evaluations, optimizer tolerance 1e-4, and TAPAS at relative gap 1e-8. Original capacities are lower bounds; upper bounds are 1.5 times those values. The illustrative investment model uses theta = 5 and budget = 100000. The objective combines equilibrium travel time and investment.

Validated measured runs: **3/3**. The designed travel time below comes from a fresh equilibrium solve at the final capacities.

| Measurement | Median |
| --- | ---: |
| Baseline total travel time | 7,480,225.310 |
| Designed total travel time | 6,186,386.519 |
| Simulated travel-time reduction | 17.30% |
| Investment used | 100,000.000 / 100,000 |
| Upper-level objective | 6,286,386.519 |
| Optimization time | 3.8049 s |
| Maximum budget violation across valid runs | 0.00e+00 |
| Maximum capacity-bound violation across valid runs | 0.00e+00 |

![Capacity optimization](figures/capacity_optimization.svg)

## Methodology and validation

- Each configuration receives one unmeasured warm-up. Full-profile TAP configurations have five measured repetitions; CNDP has three. Quick-profile configurations have one. Configuration order rotates across repetitions.
- Workers execute sequentially in fresh processes and construct fresh networks. Native computation time excludes imports, input loading, solver construction, result validation, and plotting. CNDP time includes its lower-level solves but excludes the separate baseline and fresh final-design validation.
- Statistics recording is disabled. TAPAS uses its existing seed of 42. Repeated runs measure timing variability, rather than variation across seeds.
- Independent Python checks cover finite/nonnegative flows, node conservation, BPR link costs, total travel time, and a Dijkstra-based equilibrium gap. Gap validation allows an absolute 1e-12 roundoff margin.
- CNDP also checks capacity bounds, investment accounting, and a fresh final equilibrium. Budget tolerance is 1e-6 and bound tolerance is 1e-7. Fresh and warm-started travel times must agree within 1e-5 relative tolerance.
- TAP has a 1000-iteration ceiling. Worker limits are 120 seconds for TAP and 300 seconds for CNDP, including validation. The full run has a 3600-second default deadline. Errors, timeouts, and unrun cases remain in raw results.
- Measurements follow the capacity-dependent derivative-cache correction. CNDP results describe simulated outcomes under the specified model.

## Machine and build

- CPU: 13th Gen Intel(R) Core(TM) i7-13650HX; 20 CPUs available to the process.
- Platform: Linux-7.0.0-22-generic-x86_64-with-glibc2.43.
- Compiler: x86_64-linux-gnu-g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0.
- Build: Release; release flags: -O3 -DNDEBUG.
- Python: 3.14.4; NumPy: 2.4.4.
- Git base: 6ddddf0712d451dba4ad3ea932e87f00b7d59e35; uncommitted changes present: true.
- Source fingerprint: 0431620faaf10486cef732296a632e48360c7eec23459151f1d8df73fe03c52f.
- The manifest records individual source hashes, the native extension hash, dependency versions, and exact configurations. Runtime depends on hardware, system load, compiler, and numerical stopping behavior.

## Reproduce

From the repository root, after installing the native prerequisites in the main README:

```bash
python -m pip install -e '.[examples,test]' -Ccmake.build-type=Release
python scripts/run_benchmarks.py --profile full \
  --time-budget 3600 --output performance_results/benchmark-reproduction
python scripts/run_benchmarks.py \
  --render-only performance_results/benchmark-reproduction
```

The runner downloads only the pinned TNTP inputs into build/benchmark-data. For offline preparation, add --source-root data/TransportationNetworks with a Git clone containing the pinned revision; it reads committed blobs and ignores local dataset edits. Existing run directories are protected from overwrite.

A short installation check uses --profile quick. Regenerating reports needs the plotting dependencies but does not execute solvers.

## Evidence

- [Raw measurements and detailed validation](measurements.jsonl)
- [Flat measurement table](measurements.csv)
- [Aggregated results](summary.csv)
- [Configuration and provenance](manifest.json)

Recorded unsuccessful workers (including warm-ups): **0**.
