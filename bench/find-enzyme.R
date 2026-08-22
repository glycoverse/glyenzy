devtools::load_all(quiet = TRUE)

legacy_find_enzyme_motif <- function(glycans) {
  masks <- purrr::map(
    glyenzy:::glyenzy_enzymes,
    ~ glyenzy:::.safe_have_enzyme(glycans, .x)
  )
  mask_mat <- do.call(cbind, masks)
  purrr::map(
    seq_along(glycans),
    ~ names(glyenzy:::glyenzy_enzymes)[mask_mat[.x, ]]
  )
}

glycans <- glyparse::auto_parse(c(
  "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
  "Neu5Ac(a2-8)Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
  "Gal3S(b1-3)GalNAc(a1-",
  "Gal3S6S(b1-4)GlcNAc(b1-"
))

old <- legacy_find_enzyme_motif(glycans)
new <- glyenzy:::.find_enzyme_motif(glycans)
stopifnot(identical(old, new))

time_runs <- function(f, n = 3L) {
  vapply(
    seq_len(n),
    \(i) unname(system.time(f(glycans))[["elapsed"]]),
    numeric(1)
  )
}

old_time <- time_runs(legacy_find_enzyme_motif)
new_time <- time_runs(glyenzy:::.find_enzyme_motif)
result <- data.frame(
  implementation = c("legacy", "batched"),
  median_seconds = c(median(old_time), median(new_time))
)
result$speedup <- result$median_seconds[[1]] / result$median_seconds
print(result, row.names = FALSE)
