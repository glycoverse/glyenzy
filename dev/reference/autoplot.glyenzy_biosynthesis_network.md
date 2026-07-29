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
  width = 6,
  height = 6,
  units = c("in", "cm", "mm"),
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
  names. Labels that would be too small to render after panel fitting
  are hidden with a warning.

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

- width, height:

  Positive, finite figure dimensions. Networks larger than the available
  figure are scaled proportionally, while smaller networks retain their
  natural size and are centered in the figure.

- units:

  Units for `width` and `height`. One of `"in"` (the default), `"cm"`,
  or `"mm"`.

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

A ggraph/ggplot object with a fixed-size, collision-aware layered panel
centered in the requested figure dimensions.

## Details

Glycan dimensions are measured before the network is laid out. Nodes at
the same rank are separated by their rendered bounds, while rectangular
edge caps stop arrows outside the source and target glycans. The plot
uses a fixed-size panel inside a requested figure canvas so these
clearances remain physical rather than changing with the coordinate
range.

Concrete enzyme reactions use solid edge shafts, while virtual enzyme
reactions use dashed shafts; arrowheads remain solid in both cases. This
convention is documented rather than repeated in an in-plot legend.

## Figure size

`width`, `height`, and `units` describe the complete base figure
returned by this method. Use the same dimensions with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
or the corresponding knitr/Quarto figure options. Titles, legends, or
margins added after this method returns may require additional output
space.

## Examples

``` r
if (FALSE) { # \dontrun{
network <- trace_biosynthesis_virtual(
  "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
)
network_plot <- ggplot2::autoplot(network, width = 8, height = 5)
ggplot2::ggsave(
  "biosynthesis-network.png",
  network_plot,
  width = 8,
  height = 5,
  units = "in"
)
} # }
```
