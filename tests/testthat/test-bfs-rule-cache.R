test_that("shared rules preserve selected enzyme order and duplicate labels", {
  enzymes <- lapply(
    c("B4GALT3", "MAN2A2", "B4GALT2", "MAN2A1", "B4GALT3"),
    enzyme
  )
  actual <- .bfs_search_rule_plan(enzymes)
  expected <- .prepare_bfs_rule_plan(enzymes, prepare_strict = FALSE)
  expect_identical(actual$enzyme_rule_ids, expected$enzyme_rule_ids)
  expect_equal(actual$prepared_rules, expected$prepared_rules)
  expect_equal(actual$native_rules, lapply(expected$rules, .bfs_native_rule))
})

test_that("modified rules and custom dispatch bypass shared rule preparation", {
  builtin <- enzyme("C1GALT1")
  .bfs_search_rule_plan(list(builtin))
  edited <- builtin
  edited$rules <- enzyme("B3GNT6")$rules
  custom <- builtin
  class(custom) <- c("custom_gt_enzyme", class(custom))
  for (enzyme in list(edited, custom)) {
    expect_equal(
      .bfs_search_rule_plan(list(enzyme)),
      .prepare_bfs_rule_plan(list(enzyme), prepare_strict = FALSE)
    )
  }
})
