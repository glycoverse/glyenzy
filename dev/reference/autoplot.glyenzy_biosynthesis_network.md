# Plot a glycan biosynthesis network

`autoplot.glyenzy_biosynthesis_network()` draws the typed networks
returned by
[`trace_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis.md),
[`trace_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/trace_biosynthesis_virtual.md),
[`path_biosynthesis()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis.md),
and
[`path_biosynthesis_virtual()`](https://glycoverse.github.io/glyenzy/dev/reference/path_biosynthesis_virtual.md).
A layered DAG layout determines the node positions while accounting for
converging biosynthesis routes. Parallel enzyme names are combined on
one reaction edge.

## Usage

``` r
# S3 method for class 'glyenzy_biosynthesis_network'
autoplot(
  object,
  show_enzyme = TRUE,
  size = 0.4,
  node_gap = 0.25,
  level_gap = 0.6,
  max_panel_width = 6,
  max_panel_height = 6,
  show_linkage = FALSE,
  orient = c("H", "V"),
  color_edge = FALSE,
  enzyme_label_style = c("condensed", "full"),
  ...
)
```

## Arguments

- object:

  A `glyenzy_biosynthesis_network` object returned by a biosynthesis
  function.

- show_enzyme:

  Logical. Whether to label reaction edges with enzyme or virtual-enzyme
  names.

- size:

  Positive numeric whole-cartoon scale multiplier passed to
  [`glydraw::geom_node_glycan()`](https://glycoverse.github.io/glydraw/reference/geom_node_glycan.html).
  Defaults to `0.4`.

- node_gap:

  Non-negative physical clearance, in inches, between glycan nodes at
  the same rank.

- level_gap:

  Non-negative physical clearance, in inches, between adjacent ranks.
  Increase this when enzyme labels are unusually long.

- max_panel_width, max_panel_height:

  Positive maximum panel dimensions, in inches. Networks larger than
  either limit are scaled proportionally so the complete base plot, with
  no additional outer margin, fits a graphics device of the same
  dimensions. Use `Inf` to preserve the network's natural size along one
  dimension.

- show_linkage:

  Logical. Whether glycan linkage annotations are shown. Defaults to
  `FALSE` for a compact network.

- orient:

  Glycan drawing orientation passed to
  [`glydraw::geom_node_glycan()`](https://glycoverse.github.io/glydraw/reference/geom_node_glycan.html).

- color_edge:

  Logical. Whether reaction edges and enzyme labels use the SNFG color
  of the residue added by each reaction. Defaults to `FALSE`, which
  draws them in dark grey. Concrete and virtual reactions use solid and
  dashed lines, respectively, in both modes. Reactions without one
  unambiguous added residue remain dark grey.

- enzyme_label_style:

  How parallel enzymes on one reaction edge are labeled. `"condensed"`
  (the default) groups names with the same prefix and terminal number,
  for example, `"B4GALT1/2/3, B3GALT3/4"`. `"full"` keeps complete names
  separated by `" / "`.

- ...:

  Additional glycan appearance arguments accepted by
  [`glydraw::glycanGrob()`](https://glycoverse.github.io/glydraw/reference/glycanGrob.html),
  such as `node_size`, `colors`, or `style`.

## Value

A ggraph/ggplot object with a fixed-size, collision-aware layered panel.

## Details

Glycan dimensions are measured before the network is laid out. Nodes at
the same rank are separated by their rendered bounds, while rectangular
edge caps stop arrows outside the source and target glycans. The plot
uses a fixed-size panel so these clearances remain physical rather than
changing with the coordinate range.

Concrete enzyme reactions use solid edge shafts, while virtual enzyme
reactions use dashed shafts; arrowheads remain solid in both cases. This
convention is documented rather than repeated in an in-plot legend.

## Examples

``` r
if (FALSE) { # \dontrun{
network <- trace_biosynthesis_virtual(
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
)
ggplot2::autoplot(network)
} # }
```
