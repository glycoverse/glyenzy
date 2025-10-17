# glyenzy (development version)

# glyenzy 0.2.3

## Minor improvements and bug fixes

* Updated dependency on `glymotif (>= 0.10.0)` to ensure compatibility with recent changes in package `igraph v2.2.0`.

# glyenzy 0.2.2

## Minor improvements and bug fixes

* Fix bugs introduced by the breaking changes in `glymotif` v0.7.0 and `glyrepr` v0.7.0.

# glyenzy 0.2.1

## Minor improvements and bug fixes

* Update dependencies to depend on release versions of glycoverse packages.

# glyenzy 0.2.0

## Breaking changes

* Remove `return` parameter from `find_synthesis_path()` and `rebuild_biosynthesis()`. Now they always return all possible paths.

## New features

* `rebuild_biosynthesis()` now supports multiple target glycans. The resulting graph contains all given target glycans and intermediate glycans.

## Minor improvements and fixes

* Update documentations of `find_synthesis_path()` and `rebuild_biosynthesis()` to include some important notes.
* Fix bugs in `find_synthesis_path()` and `rebuild_biosynthesis()` that glycan structure strings other than IUPAC-condensed format cannot be parsed.
* Add checks in all functions to ensure that the input glycans have intact linkages and no substituents, and raise errors with helpful messages if not.
* Refactor `apply_enzyme()` to stop using internal `glymotif` functions to avoid fragile dependency.

# glyenzy 0.1.1

## Minor improvements and fixes

* Fix bugs introduced by the breaking changes in `glyrepr` v0.7.0.

# glyenzy 0.1.0

First release of `glyenzy`!
