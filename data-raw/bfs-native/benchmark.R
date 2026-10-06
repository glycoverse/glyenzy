source("data-raw/bfs-native/helpers.R")
source("data-raw/bfs-native/cases.R")
selection <- Sys.getenv("BFS_CASE", "")
if (nzchar(selection)) {
  cases <- cases[selection]
}
results <- list()
checks <- list()
for (label in names(cases)) {
  cat("Case", label, "\n")
  args <- cases[[label]]
  w <- do.call(prepare, args)
  expected <- normalize_ref(reference(w))
  actual <- decode(do.call(native_bfs, w$native), w)
  checks[[label]] <- parity(expected, actual)
  print(checks[[label]])
  if (!all(checks[[label]])) {
    stop("Semantic or ordering parity failed")
  }
  for (i in seq_len(7L)) {
    for (arm in if (i %% 2L) c("R", "native") else c("native", "R")) {
      gc()
      timing <- system.time({
        if (arm == "R") {
          out <- normalize_ref(reference(w))
        } else {
          fresh <- do.call(prepare, args)
          out <- decode(do.call(native_bfs, fresh$native), fresh)
        }
      })[["elapsed"]]
      stopifnot(all(parity(expected, out)))
      results[[length(results) + 1L]] <- data.frame(
        case = label,
        run = i,
        arm = arm,
        elapsed = timing,
        edges = nrow(out$edges),
        visited = length(out$parent) + 1L
      )
      cat(arm, i, timing, "\n")
    }
  }
}
write.csv(
  do.call(rbind, results),
  file.path(base_dir, "timings.csv"),
  row.names = FALSE
)
saveRDS(checks, file.path(base_dir, "parity.rds"))
capture.output(sessionInfo(), file = file.path(base_dir, "session.txt"))
