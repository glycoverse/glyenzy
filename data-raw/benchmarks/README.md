# Biosynthesis tracing benchmark

Tested on 2026-10-19 (Asia/Shanghai), using macOS Tahoe 26.5.2, ARM64, and R 4.6.1.

Baseline: `e925d5f` (`glyenzy` 0.10.0.9000). Optimized: the implementation committed with this report. Both were installed into separate temporary libraries with the normal `-O2` build settings and the same dependencies.

## Results

| Case | Targets | Baseline mean (ms) | Optimized mean (ms) | Optimized median (ms) | Optimized p95 (ms) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| high_mannose_N | 1 | 459.06 | 29.64 | 29.00 | 31.55 | 15.49x |
| complex_N | 1 | 622.22 | 36.52 | 37.00 | 38.00 | 17.04x |
| O_GalNAc | 1 | 355.12 | 24.66 | 24.00 | 29.75 | 14.40x |
| O_Man | 1 | 354.22 | 23.76 | 23.00 | 26.00 | 14.91x |
| example_single | 1 | 462.62 | 24.98 | 24.50 | 28.00 | 18.52x |
| example_multi | 2 | 699.64 | 28.90 | 29.00 | 31.00 | 24.21x |

The five single-target cases have an equal-weight mean of **27.912 ms**, compared with **450.648 ms** for the baseline (16.15x faster). Every single-target case has a mean below 50 ms. The two-target case is timed per complete call, not per target.

## Measurement method

- Default arguments and character-string input, including parsing, enzyme selection, full search, and final network construction.
- One initial call per case, followed by 50 measured calls per case in randomly interleaved order (fixed seed). Natural garbage collection is included. Package loading and compilation are excluded.
- Baseline and optimized runs were sequential, after tests and compilation finished.
- These are warmed-session measurements. The first call in the optimized session took **347 ms** (baseline: 514 ms), including one-time preparation of the built-in rules. Later cases share that preparation; their `first_ms` values are not independent cold starts.
- Only target-independent rule preparation is cached. Targets, search states, and final networks are recomputed on every call.
- Results describe these exact inputs on this machine; glycan complexity affects search cost.

## Validation

- All 3,397 test assertions passed, with no failures, warnings, or skips.
- `R CMD check --no-manual --no-tests` completed with 0 errors and 0 warnings; tests had already run separately. Its two environment notes concern unavailable remote clock verification and the macOS `xcrun_db` temporary cache. Examples and vignette builds passed.
- Added differential tests for concrete/generic and intact/topological inputs, sulfate subsets, N-glycan trimming, custom enzymes, extra graph metadata, and fallback behavior.
- Baseline and optimized installed builds returned identical classes, graph attributes, ordered vertices, ordered edges, and enzyme annotations for all six benchmark cases.

## Reproduce

Install each revision into its own library, then run sequentially:

```sh
Rscript data-raw/benchmarks/trace-biosynthesis.R /path/to/baseline-library /tmp/trace-baseline 50
Rscript data-raw/benchmarks/trace-biosynthesis.R /path/to/optimized-library /tmp/trace-optimized 50
```

The script contains the exact input structures and writes summary statistics, every timed call, session information, and network snapshots for equality checks. This directory retains the measured summaries, all 600 timings, and the optimized session information.
