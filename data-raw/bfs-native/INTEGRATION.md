# Production integration

The package now uses `cpp_bfs_search()` for every BFS search. The R6 wrapper
validates/prepares inputs, supplies R extension callbacks, and restores result
objects. The old R expansion, integration and queue loops have been removed
from package code. Their frozen implementation exists only in the test fixture
`tests/testthat/fixtures/bfs-reference.R` for differential testing.

## Execution paths

| Path | Execution |
|---|---|
| Ordinary and abstract GT/GH | Compact native matching, actions, pruning, canonicalization, cache and queue |
| ST and sequential sulfation | Native sulfate/occupied-position checks and substituent matching |
| Concrete/generic, intact/topological, partial internal inputs | Native strict/lenient matching and topology reduction; public input validation unchanged |
| Key and whole target matching | Native target matching, with duplicate hits and per-level completion retained |
| N-precursor pruning | Native pre-MGAT2 detection and terminal Glc/Man/Hex trimming |
| Enzyme glycan-type restrictions | Native classification for each substrate, including changes during traversal |
| Starter/precursor-transferase entries | Same inactive BFS behavior as the reference |
| Custom S3 enzymes | Original R dispatch per substrate/enzyme cell; native traversal/bookkeeping |
| Stateful filters | Original R filter and logical-vector validation per nonempty, pruned cell; native traversal/bookkeeping |
| Attribute-bearing, floating or noncompact graphs/rules | Original R graph operations through explicit callbacks; native traversal/bookkeeping |
| Virtual-step fallback | Existing orchestration calls the same native BFS entry point |

There is no switch back to an R BFS implementation. Retaining R callbacks for
user R functions and specialized graph operations is deliberate compatibility,
not an unimplemented search path. Canonical sorting uses `base::order` outside
C collation where required; user locale is not silently replaced by byte order.
Graph imports validate rooted-tree connectivity before native recursive work.

## Verification

Verified implementation commit: `f87c263214cc9e5545cd19e9adf2ad96d314547e`.
The final macOS R 4.6.1 `R CMD check --as-cran --no-manual` completed with
**0 errors, 0 warnings, 0 notes**. The full suite passed **3,319 assertions**,
with no failures, warnings or skips. Examples and vignette rebuilding passed.
Raw evidence is retained in integrated-check.log, integrated-tests.log and
integrated-install.log. The final optimized binary and source checksums are
recorded in integrated-build.txt.

- Differential tests compare ordered reaction edges, target hits (including
  duplicates/order), parent maps, parent enzymes, steps and missing target keys.
- Every bundled ordinary/abstract enzyme rule is compared in intact and
  topological mode against the frozen R engine.
- Additional cases cover multistep GT/GH, multiple N targets, sulfate subsets,
  generic residues, filter call order/validation, locale, alditol, graph metadata,
  floating metadata, initial targets, unreachable targets, repeated goals,
  empty enzymes, zero/fractional/infinite internal step bounds and partial search.
- Existing public API tests exercise custom S3 dispatch, mixed enzyme cells,
  public input rejection and virtual-step recovery.

## Installed-build benchmark

Run from the repository root after installing an optimized build in an isolated
library:

```sh
R CMD INSTALL --preclean --library=/tmp/glyenzy-native-lib .
BFS_NATIVE_LIB=/tmp/glyenzy-native-lib Rscript data-raw/bfs-native/benchmark-integrated.R
```

The library directory must exist. Results are stored separately from the
historical standalone prototype in `integrated-timings.csv`,
`integrated-summary.csv` and `integrated-session.txt`.

Both arms start with the same parsed inputs and enzymes. Timing includes the
complete internal BFS call: setup, conversion, search and result restoration.
Compilation, input parsing and public network assembly are excluded. Seven
repetitions alternate arm order; each result must exactly match the frozen R
engine before its timing is retained. This measures the implemented native BFS,
not the entire `trace_biosynthesis()` API and not a language-only effect.

### Final timings

| Case | Frozen R median | Native median | Speedup |
|---|---:|---:|---:|
| Branched O-glycan | 150 ms | 18 ms | 8.33x |
| Complex N-glycan | 180 ms | 16 ms | 11.25x |
| N-precursor GT/GH | 429 ms | 38 ms | 11.29x |
| Multi-target N-glycans | 4.134 s | 65 ms | 63.60x |

All seven paired runs favored native for every case, with exact ordered-output
and parent-map equality in every repetition. These are four explicit workloads,
not a guarantee of the same speedup for arbitrary workloads or R callbacks.
The integrated results are not directly comparable to the earlier prototype
ratios: that experiment charged the prototype extra input parsing/preparation
and rebuilt every visited graph; this comparison measures the actual production
entry point with equal prepared inputs.
