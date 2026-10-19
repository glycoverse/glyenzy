# Bundled native rules shared by all target sets. Compare full enzyme objects
# so edited rules and custom S3 classes retain runtime preparation. The stored
# records contain plain R values, with no serialized igraph instances.
.bfs_builtin_rule_plan <- function() {
  glyenzy_bfs_rules
}

.build_bfs_builtin_rule_plan <- function(enzymes) {
  plan <- .prepare_bfs_rule_plan(enzymes, prepare_strict = FALSE)
  list(
    enzymes = lapply(enzymes, function(enzyme) {
      structure(
        list(name = enzyme$name, glycan_type = enzyme$glycan_type),
        class = c("glyenzy_bfs_enzyme", class(enzyme))
      )
    }),
    default_indices = which(
      !vapply(
        enzymes,
        \(enzyme) .is_starter_gt(enzyme) || .is_npre_gt(enzyme),
        logical(1)
      )
    ),
    enzyme_rule_ids = plan$enzyme_rule_ids,
    native_rules = lapply(plan$rules, .bfs_native_rule)
  )
}

.bfs_search_rule_plan <- function(enzymes) {
  builtin <- vapply(
    enzymes,
    function(enzyme) {
      .bfs_is_builtin_enzyme(enzyme)
    },
    logical(1)
  )
  if (!all(builtin) || !length(enzymes)) {
    return(.prepare_bfs_rule_plan(enzymes, prepare_strict = FALSE))
  }
  plan <- .bfs_builtin_rule_plan()
  indices <- match(
    vapply(enzymes, `[[`, character(1), "name"),
    names(plan$enzymes)
  )
  rule_ids <- plan$enzyme_rule_ids[indices]
  selected <- unique(unlist(rule_ids, use.names = FALSE))
  list(
    enzyme_rule_ids = lapply(rule_ids, match, table = selected),
    native_rules = plan$native_rules[selected]
  )
}

.bfs_is_builtin_enzyme <- function(enzyme) {
  if (is.null(enzyme$name)) {
    return(FALSE)
  }
  if (inherits(enzyme, "glyenzy_bfs_enzyme")) {
    return(identical(enzyme, glyenzy_bfs_rules$enzymes[[enzyme$name]]))
  }
  identical(enzyme, glyenzy_enzymes[[enzyme$name]])
}

.bfs_resolve_enzyme <- function(enzyme) {
  if (inherits(enzyme, "glyenzy_bfs_enzyme")) {
    return(glyenzy_enzymes[[enzyme$name]])
  }
  enzyme
}

.bfs_default_enzymes <- function(targets) {
  records <- lapply(
    glyrepr::get_structure_graphs(targets, return_list = TRUE),
    .bfs_native_graph
  )
  if (any(vapply(records, is.null, logical(1)))) {
    return(NULL)
  }
  glyenzy_bfs_rules$enzymes[glyenzy_bfs_rules$default_indices]
}
