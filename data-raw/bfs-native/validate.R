source("data-raw/bfs-native/helpers.R")
source("data-raw/bfs-native/cases.R")
# Empty-edge cases still compare typed, ordered result tables.
original_normalize <- normalize_ref
normalize_ref <- function(z) {
  x <- original_normalize(z)
  if (is.null(x$edges)) {
    x$edges <- data.frame(
      from = character(),
      to = character(),
      enzyme = character(),
      step = integer()
    )
  }
  for (n in c("parent", "parent_enzyme")) {
    if (is.null(x[[n]])) x[[n]] <- setNames(character(), character())
  }
  if (is.null(x$parent_step)) {
    x$parent_step <- setNames(integer(), character())
  }
  x
}
cases$zero_steps <- list("GalNAc(a1-", "Gal(b1-3)GalNAc(a1-", "C1GALT1", 0L)
cases$initial_target <- list("GalNAc(a1-", "GalNAc(a1-", "C1GALT1", 3L)
cases$unreachable <- list("GalNAc(a1-", "GlcNAc(b1-3)GalNAc(a1-", "C1GALT1", 3L)
cases$duplicate_targets <- list(
  "GalNAc(a1-",
  rep("Gal(b1-3)GalNAc(a1-", 2),
  c("C1GALT1", "C1GALT1"),
  1L
)
cases$partial <- cases$n_complex
cases$partial[[4]] <- 2L
checks <- list()
outcomes <- list()
for (label in names(cases)) {
  w <- do.call(prepare, cases[[label]])
  a <- normalize_ref(reference(w))
  b <- decode(do.call(native_bfs, w$native), w)
  checks[[label]] <- parity(a, b)
  outcomes[[label]] <- data.frame(
    case = label,
    targets = length(unique(as.character(w$to))),
    reached = length(unique(b$found)),
    missing = b$missing,
    edges = nrow(b$edges),
    visited = length(b$parent) + 1L
  )
  cat(
    label,
    paste(names(checks[[label]])[!checks[[label]]], collapse = ","),
    "\n"
  )
  stopifnot(all(checks[[label]]))
}
saveRDS(checks, "data-raw/bfs-native/validation.rds")

write.csv(
  do.call(rbind, outcomes),
  "data-raw/bfs-native/validation_outcomes.csv",
  row.names = FALSE
)
