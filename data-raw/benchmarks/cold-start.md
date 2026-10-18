# First-call biosynthesis tracing benchmark

Measured on 2026-10-19 (Asia/Shanghai), macOS Tahoe 26.5.2, ARM64, R 4.6.1.
Baseline: `33a70cd`, the previous warm-session optimization. Optimized: the
implementation committed with this report. Both are normal installed builds
with the same dependency libraries and compiler settings.

## First-call results

Each sample starts a fresh `Rscript --vanilla` process, loads `glyenzy`, then
calls `trace_biosynthesis()` exactly once with default arguments and a character
input. There is no pre-GC, warm-up call, or target/result cache. Parsing, enzyme
selection, full search, final network construction, lazy initialization, and any
natural garbage collection inside the call are included. Timing uses elapsed
`proc.time()` with approximately 1 ms resolution. Processes run sequentially;
case and revision order are randomized with a fixed seed. Compilation, tests,
and other benchmark runs finished before these measurements.

| Case | Fresh processes per version | Previous mean (ms) | New mean (ms) | New median (ms) | New p95 (ms) | New maximum (ms) |
|---|---:|---:|---:|---:|---:|---:|
| high_mannose_N | 20 | 343.85 | 30.35 | 30.0 | 32.05 | 33 |
| complex_N | 20 | 345.65 | 36.70 | 37.0 | 38.00 | 38 |
| O_GalNAc | 20 | 334.35 | 28.50 | 28.5 | 30.00 | 30 |
| O_Man | 20 | 333.80 | 28.95 | 29.0 | 30.00 | 30 |
| example_single | 20 | 333.75 | 29.45 | 30.0 | 31.05 | 32 |

The equal-weight single-target first-call mean fell from **338.28 ms to 30.79 ms**
(10.99x faster). All 100 optimized first calls were between **27 and 38 ms**.
All 200 networks, including both versions and every repetition, were identical
in classes, graph attributes, ordered vertices/edges, and enzyme annotations.
These figures describe the five exact input structures in the script on this
machine; search cost still depends on glycan complexity and chosen options.

Package loading is measured separately and is **excluded** from the function
timing above. Mean `library(glyenzy)` time changed from **329.35 to 355.53 ms**;
imports now load the namespaces used by every default search. No search is run
at package load. Mean package load plus first call improved from **667.63 to
386.32 ms**. R process startup is excluded from both figures.

## Warm-session regression

The existing benchmark script was also run for 50 randomly interleaved repetitions
per case, including natural garbage collection. The five single-target means
average **15.168 ms**. All six networks, including the two-target example,
exactly match the snapshots from `33a70cd`.

| Case | Targets | Mean (ms) | p95 (ms) |
|---|---:|---:|---:|
| high_mannose_N | 1 | 14.42 | 16.00 |
| complex_N | 1 | 21.00 | 22.55 |
| O_GalNAc | 1 | 13.00 | 15.00 |
| O_Man | 1 | 13.72 | 15.00 |
| example_single | 1 | 13.70 | 16.55 |
| example_multi | 2 | 15.40 | 17.00 |

## Implementation

- Generate native rule records, motif plans, and fixed starting structures with
  the enzyme database. The database-generation script remains the source for
  `R/sysdata.rda`; regression tests check that generated plans match the database.
- Default native searches use compact enzyme descriptors. Full enzyme objects
  and R rule preparation are loaded on demand for graph-metadata callbacks.
  Custom enzymes, explicit filters, and virtual fallback retain their existing
  processing paths.
- Parse canonical condensed strings directly, retaining automatic format
  detection and diagnostics for other inputs. Use a lightweight intact-structure
  check, falling back to the public classifier for ambiguous or floating input.
- Construct result tables without loading the tibble construction helpers on
  the first call. Each target and search result is recomputed.

## Validation

- Full test suite: **3,430 passing assertions**, zero failures, warnings, or skips.
- `R CMD check --no-manual --no-tests`: zero errors and warnings; tests ran
  separately above. Examples and vignette builds passed. The two environment
  notes concern unavailable remote clock verification and the macOS `xcrun_db`
  temporary cache.
- Added coverage for generated-plan consistency, default tracing without full
  enzyme graphs, fallback callbacks with metadata, direct-parser parity, and
  intact/ambiguous/topological/missing/floating structure classification.
- Installed-build comparisons checked all five cold-start cases in every
  process, plus the sixth, multi-target case in the warm benchmark.

## Reproduce

Install `33a70cd` and the optimized revision into separate libraries, then run:

```sh
Rscript data-raw/benchmarks/trace-biosynthesis-cold.R /path/to/previous-library /path/to/optimized-library /tmp/trace-cold 20
Rscript data-raw/benchmarks/trace-biosynthesis.R /path/to/optimized-library /tmp/trace-warm 50
```

The cold script contains the exact inputs and fails if any network differs
between revisions or repetitions. This directory retains the 200 fresh-process
measurements, per-case summaries, and session information for both versions.
Warm-session measurements and summaries from the follow-up are retained
separately from the [original warm benchmark](README.md).
