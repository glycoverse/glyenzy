# Changelog

## glyenzy (development version)

### New features

- Add first-class sulfotransferases (`ST`) with 12 human N- and O-glycan
  enzymes, sulfate-aware inference and biosynthesis, and virtual
  sulfation steps. (#20, \#33)
- [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  and
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  gain `max_virtual_steps` to bridge a bounded number of unsupported,
  target-directed transitions before resuming concrete enzyme tracing;
  fallback edges are marked by `is_virtual`. (#32)
- New
  [`trace_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis_virtual.md)
  and
  [`path_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis_virtual.md)
  build enzyme-agnostic networks by trimming targets backward;
  `annotate_enzymes = TRUE` adds exact rule-matched candidates in
  `concrete_enzymes`. (#22, \#28, \#30)
- Add support for non-intact glycan structures by using lenient motif
  matching with a warning about reduced reliability. (#25)

### Minor improvements and fixes

- Speed up large and multi-target biosynthesis searches by vectorizing
  enzyme prefiltering, sharing frontier motif matching, caching
  equivalent graph products, pruning irreversible pre-MGAT2 decorations,
  prioritizing inclusive targets, and using one multi-target
  reachability traversal. (#34)
- Biosynthesis searches now keep intermediate products as graphs through
  pruning and only generate canonical structure keys for candidates that
  enter the network, using the low-level graph APIs in `glyrepr` 0.14.0.
  (#29)
- [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  now reuses prepared glycan graphs, shares equivalent enzyme-rule work,
  batches breadth-first expansion, and rejects products that would reuse
  an occupied acceptor carbon. (#26)
- To ensure robust biosynthesis network inference, paucimannose
  N-glycans are not supported anymore. (#25)

## glyenzy 0.6.3

### Minor improvements and fixes

- [`apply_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/apply_enzyme.md),
  [`have_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/have_enzyme.md),
  [`count_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/count_enzyme.md),
  and biosynthesis-path helpers now work with `glymotif` 0.17.0 and
  later.

## glyenzy 0.6.2

### Minor improvements and fixes

- Update built-in enzyme data and require `glyrepr` 0.13.0 or later for
  compatibility with refreshed glycan structure data (#24).

## glyenzy 0.6.1

### Minor improvements and fixes

- Fix rules for MAN1A1, MAN1A2, and MAN1C1.

## glyenzy 0.6.0

### New features

- Add `method = "path"` to
  [`have_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/have_enzyme.md),
  [`count_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/count_enzyme.md),
  [`match_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/match_enzyme.md),
  and
  [`find_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/find_enzyme.md)
  for biosynthesis-derived enzyme inference. (#19)

### Minor improvements and fixes

- Fix enzyme printing without requiring users to attach `glyrepr`.
  (f0e42b2)

## glyenzy 0.5.4

### Minor improvements and fixes

- Fix
  [`match_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/match_enzyme.md)
  to ignore enzyme reject motifs when reporting matched residues,
  consistent with
  [`have_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/have_enzyme.md)
  and
  [`count_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/count_enzyme.md).
  (0f18c74)

## glyenzy 0.5.3

### Minor improvements and fixes

- Update rules for FUT3, FUT4, FUT5, FUT6, FUT7, and FUT9. (2c6b28b)
- Remove rules for B4GALT5, B4GALT6, and ST3GAL5 because these enzymes
  are exclusive for glycolipids. (798ea6f, 8dd30c8)
- Update ST3GAL rules to reject Lewis antigen acceptors. (c969912)

## glyenzy 0.5.2

### Minor improvements and fixes

- Fix rules for FUT3, FUT4, FUT5, FUT6, and FUT9. (5e7d657)

## glyenzy 0.5.1

### Minor improvements and fixes

- Fix rules for FUT1, FUT2, and B3GALT5. (50f3d5d)

## glyenzy 0.5.0

### Breaking changes

- Rename APIs to use shorter, consistent names (#2):
  - `all_enzymes()` to
    [`db_enzymes()`](https://glycoverse.github.io/glyenzy/dev/reference/db_enzymes.md)
  - `create_enzyme()` to
    [`make_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/make_enzyme.md)
  - `get_involved_enzymes()` to
    [`find_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/find_enzyme.md)
  - `is_synthesized_by()` to
    [`have_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/have_enzyme.md)
  - `count_enzyme_steps()` to
    [`count_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/count_enzyme.md)
  - `find_synthesis_path()` to
    [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  - `rebuild_biosynthesis()` to
    [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  - `spawn_glycans()` to
    [`grow_glycans()`](https://glycoverse.github.io/glyenzy/dev/reference/grow_glycans_step.md)
  - `spawn_glycans_step()` to
    [`grow_glycans_step()`](https://glycoverse.github.io/glyenzy/dev/reference/grow_glycans_step.md)

### New features

- Add support for initiating enzymes including DPAGT1, FUT10, FUT11,
  POFUT1, POFUT2, POGLUT1, POGLUT2, POGLUT3, POMT1, POMT2, TMTC1, TMTC2,
  TMTC3, TMTC4, and GALNT1 through GALNT19. (#5)
- Add support for N-glycan precursor synthesis enzymes (ALGs). (#6)
- Add
  [`match_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/match_enzyme.md)
  to identify the residues added by a glycosyltransferase. (#3)
- Add
  [`view_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/view_enzyme.md)
  to visualize residues added by an enzyme on a glycan cartoon. (#4)

### Minor improvements and fixes

- Fix the bug that find_enzyme didn’t support paucimannose glycans.
  (1c79a04)

## glyenzy 0.4.3

### Minor improvements and fixes

- Update dependency strategy to use the r-universe repo.

## glyenzy 0.4.2

### Minor improvements and fixes

- Update internal data to adapt to glyrepr 0.10.0.

## glyenzy 0.4.1

### Minor improvements and fixes

- glyenzy now depends on the CRAN version of glyparse.

## glyenzy 0.4.0

### Breaking changes

- Strengthen validation of the `rejects` field in enzyme rules: each
  reject must contain the acceptor motif.
- Rework `rejects` handling in
  [`apply_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/apply_enzyme.md)
  (and callers) to evaluate rejects per match rather than per glycan,
  improving accuracy for glycans with multiple acceptor matches.
- Update built-in rules for MAN2A1, MAN2A2, MGAT3, MGAT4A, and MGAT4B to
  align with the new reject semantics.
- Remove the now-redundant `rejects_alignment` field from `enzyme_rule`
  objects; `acceptor_alignment` is reused for reject checks.

### New features

- Add
  [`make_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/make_enzyme.md),
  enabling programmatic construction of custom enzymes.

### Minor improvements and fixes

- [`grow_glycans()`](https://glycoverse.github.io/glyenzy/dev/reference/grow_glycans_step.md)
  now shows live glycan counts in the progress bar.

## glyenzy 0.3.2

### Minor improvements and fixes

- Fix errors in the rules of some enzymes, including B3GALT1, B3GALT2,
  and FUT8.
- Fix a bug in progress bar of
  [`grow_glycans()`](https://glycoverse.github.io/glyenzy/dev/reference/grow_glycans_step.md).

## glyenzy 0.3.1

### Minor improvements and fixes

- glyenzy now depends on the CRAN version of glyrepr.

## glyenzy 0.3.0

### Breaking changes

- Remove FUT10, as it belongs to the starting point of some O-Fuc
  glycans.

### New features

- Add many new enzymes, including FUT5, FUT6, ST3GAL5, ST6GALNAC5,
  ST6GALNAC6, ST8SIA1, ST8SIA5, ST8SIA6, B3GAT3, CHPF, CHPF2, CHSY1,
  CHSY3, EXT1, EXT2, HAS1, HAS2, HAS3, LARGE1, LARGE2, B3GLCT, A4GALT,
  ABO, B3GALT6, B4GALT6, GXYLT1, GXYLT2, LARGE1, LARGE2, XXYLT1,
  B3GALNT1, B4GALNT1, B4GALNT3, B4GALNT4, CSGALNACT1, CSGALNACT2,
  B3GNT5, B3GNT7, POMGNT1, POMGNT2, LFNG, MFNG, RFNG, EXLT1, EXLT2,
  EXLT3
- Update the rules of some existing enzymes, including FUT1, FUT2, FUT3,
  FUT4, FUT7, FUT8, FUT9, ST3GAL2, ST3GAL3, ST3GAL4, ST6GALNAC1,
  ST6GALNAC3, ST8SIA2, ST8SIA3, ST8SIA4, B3GALT4, B3GALT5, B4GALT1,
  B4GALT4, B4GALT5, B4GALNT2, MGAT5B, B3GNT8, GCNT2, GCNT3, A4GNT,
  ST3GAL1
- [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  now supports many new glycan types, including O-Man, O-GlcNAc, O-Fuc,
  and O-Glc.

### Minor improvements and fixes

- [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  and
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  has been optimized for large glycans.
- All functions explicitly check if the input glycans are concrete
  (e.g. “Glc”, “GalNAc”) and raise errors with helpful messages if not.

## glyenzy 0.2.3

### Minor improvements and bug fixes

- Updated dependency on `glymotif (>= 0.10.0)` to ensure compatibility
  with recent changes in package `igraph v2.2.0`.

## glyenzy 0.2.2

### Minor improvements and bug fixes

- Fix bugs introduced by the breaking changes in `glymotif` v0.7.0 and
  `glyrepr` v0.7.0.

## glyenzy 0.2.1

### Minor improvements and bug fixes

- Update dependencies to depend on release versions of glycoverse
  packages.

## glyenzy 0.2.0

### Breaking changes

- Remove `return` parameter from
  [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  and
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md).
  Now they always return all possible paths.

### New features

- [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  now supports multiple target glycans. The resulting graph contains all
  given target glycans and intermediate glycans.

### Minor improvements and fixes

- Update documentations of
  [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  and
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  to include some important notes.
- Fix bugs in
  [`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md)
  and
  [`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md)
  that glycan structure strings other than IUPAC-condensed format cannot
  be parsed.
- Add checks in all functions to ensure that the input glycans have
  intact linkages and no substituents, and raise errors with helpful
  messages if not.
- Refactor
  [`apply_enzyme()`](https://glycoverse.github.io/glyenzy/dev/reference/apply_enzyme.md)
  to stop using internal `glymotif` functions to avoid fragile
  dependency.

## glyenzy 0.1.1

### Minor improvements and fixes

- Fix bugs introduced by the breaking changes in `glyrepr` v0.7.0.

## glyenzy 0.1.0

First release of `glyenzy`!
