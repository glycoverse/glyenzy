# Fused native BFS experiment (2026-10-06)

This is an isolated feasibility experiment. Package code, public APIs,
dependencies, and the default BFS backend are unchanged.

## Results

Seven alternating repetitions per arm, elapsed seconds (median):

| Workload | Edges / visited | R BFS | Native with conversion | Speedup |
|---|---:|---:|---:|---:|
| Branched O-glycan | 16 / 9 | 0.173 | 0.045 | 3.84x |
| Complex N-glycan | 23 / 13 | 0.193 | 0.034 | 5.68x |
| N-precursor, GT + GH | 37 / 18 | 0.443 | 0.061 | 7.26x |
| Multi-target N-glycans | 344 / 93 | 4.427 | 0.180 | 24.59x |

All seven paired runs favored native in each case. All eight parity assertions
passed for all four workloads and five boundary cases. The machine was not
isolated: the multi-target timings varied from 4.067–9.657 s for R and
0.165–0.515 s for native; use medians and retain these ranges when interpreting
precision. The speedup is substantial on these workloads but is not a guarantee
for every enzyme set, glycan, or the public tracing API.

Conclusion: a fused native graph representation is worth pursuing. This experiment
combines a C++ implementation with fewer graph conversions, deferred canonical
string construction, and compact structural keys. It does not isolate a language-
only effect. Complete compatibility and fallback routing remain production work.

## Reproduce

From the glyenzy repository root, with Rcpp and BH installed:

```sh
Rscript -e 'devtools::load_all(); source("data-raw/bfs-native/benchmark.R"); source("data-raw/bfs-native/validate.R"); source("data-raw/bfs-native/summarize.R")'
```

`BFS_CASE=n_multi` optionally selects one timing workload. Running a subset
replaces timings.csv with that subset; omit it to reproduce the full report.
Compilation and warm-up are outside the measured intervals.

## What runs in C++

One call performs rule matching, requirements and site-specific rejects,
GT addition / terminal GH removal, target containment and pre-MGAT2 pruning,
structural caching, per-enzyme product deduplication, visited-set lookup,
parent and reaction-edge recording, target detection, and queue progression.
Products stay in compact native graph profiles; igraph objects and canonical
IUPAC strings are built only when returning results. There are no R function
callbacks during native search. The matcher still allocates Rcpp matrices/lists
and uses the same Boost VF2 algorithm as glymotif; this is not an allocation-free
or parallel implementation.

Canonical depth/linkage/subtree-signature traversal is essential: an earlier
prototype preserved the reaction set but changed edge ordering and first parents
in the multi-target case. The final version normalizes newly queued nodes and
requires exact ordering and parent parity before timing.

## Comparison and limits

Reference: current `bfs_synthesis_search()` at glyenzy commit
`65b346871dee612cdf044a16b38d7f449c860d98` (already batched and cached).
See session.txt for R, platform, locale, and dependency versions; see
vendor/PROVENANCE.md for the frozen native matcher source.

Each case uses seven repetitions per arm, alternating R/native and native/R.
Garbage collection occurs before each arm. The R arm includes engine setup,
search, and extraction of result tables from its environments. The native arm
additionally includes reparsing inputs, loading enzymes, preparing compact
rules/targets, the search, rebuilding output graphs, canonical IUPAC strings,
and the equivalent result tables. Thus the native timings include conversion
costs and conservatively charge it extra input preparation; these are not pure
C++ kernel timings. Public `trace_biosynthesis()` prefiltering and final network
assembly are not measured. No claim of the same speedup for that entire API is
made.

Before timing and on every repetition, the harness checks exact reaction edges
including their order, reached targets including multiplicity/order, first
parents, parent enzymes, discovery steps, and missing-target count. Boundary
validation additionally covers zero steps, source already at target, an
unreachable target, duplicate targets/enzymes, and a truncated search.

Supported experimental scope is intact concrete, unsubstituted, nonfloating,
non-alditol glycans with standard GT/GH rules and key-based target matching.
The wrapper rejects enzyme glycan-type restrictions incompatible with the
source. It assumes that compatible types remain compatible during these
workloads. Sulfotransferases, topological/partial structures, floating structures,
custom S3 enzymes, stateful filters, whole/lenient target matching, warnings/error
API compatibility, and cross-locale canonical tie ordering are not validated.
This is not a drop-in production backend.

The four workload inputs are explicit in cases.R. They test O-glycan branching,
complex N-glycan elongation, N-precursor trimming with GH/GT, and multi-target
N-glycan branching with competing enzymes. They are selected mechanistic cases,
not a population sample of all glycans or a replay of the historical 55-target
workload. Reported speedups are practical measured ratios, not a claim of
statistical generalization.

Results: summary.csv (medians, ranges, speedup, graph sizes), timings.csv (raw
per-run seconds), parity.csv and validation.csv (individual assertions).
