The matching implementation is copied from glymotif commit
`f1d7d3a718435f70e743e1c47e154cd39c1a96ad`:

- `vf2.h`: complete `src/vf2.cpp`.
- `matcher.h`: lines 1–173 of `src/structure-matcher.cpp`, containing the compact
  graph profile, residue compatibility and matching implementation.

MIT license: see LICENSE.md and LICENSE. These files are frozen experimental
inputs, not a new public interface or dependency on an installed C++ ABI.

The canonical traversal in native.cpp mirrors the depth/linkage/subtree-signature
ordering in glyrepr's `R/structure-to-iupac.R`. It only handles the explicitly
guarded unsubstituted, nonfloating trees used by this experiment. The C locale's
collation is used for signature ties; cross-locale equivalence is not established.
