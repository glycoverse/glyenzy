# A fixed-size cache of the bundled rules, shared by all target sets. Search
# states and results are never cached. Compare full enzyme objects so edited
# rules and custom S3 classes always take the uncached preparation path.
.bfs_builtin_rule_plan <- local({
  source <- NULL
  plan <- NULL
  function() {
    if (is.null(plan) || !identical(source, glyenzy_enzymes)) {
      plan <<- .prepare_bfs_rule_plan(glyenzy_enzymes, prepare_strict = FALSE)
      plan$native_rules <<- lapply(plan$rules, .bfs_native_rule)
      source <<- glyenzy_enzymes
    }
    plan
  }
})

.bfs_search_rule_plan <- function(enzymes) {
  builtin <- vapply(
    enzymes,
    function(enzyme) {
      !is.null(enzyme$name) && identical(enzyme, glyenzy_enzymes[[enzyme$name]])
    },
    logical(1)
  )
  if (!all(builtin) || !length(enzymes)) {
    return(.prepare_bfs_rule_plan(enzymes, prepare_strict = FALSE))
  }
  plan <- .bfs_builtin_rule_plan()
  indices <- match(
    vapply(enzymes, `[[`, character(1), "name"),
    names(glyenzy_enzymes)
  )
  rule_ids <- plan$enzyme_rule_ids[indices]
  selected <- unique(unlist(rule_ids, use.names = FALSE))
  list(
    rules = plan$rules[selected],
    prepared_rules = stats::setNames(
      plan$prepared_rules[indices],
      names(enzymes)
    ),
    enzyme_rule_ids = lapply(rule_ids, match, table = selected),
    native_rules = plan$native_rules[selected]
  )
}
