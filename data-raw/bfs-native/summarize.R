x <- read.csv("data-raw/bfs-native/timings.csv")
summary <- do.call(
  rbind,
  lapply(split(x, x$case), function(z) {
    r <- z$elapsed[z$arm == "R"]
    cpp <- z$elapsed[z$arm == "native"]
    data.frame(
      case = z$case[[1]],
      edges = z$edges[[1]],
      visited = z$visited[[1]],
      R_median = median(r),
      native_median = median(cpp),
      speedup = median(r) / median(cpp),
      R_min = min(r),
      R_max = max(r),
      native_min = min(cpp),
      native_max = max(cpp),
      paired_wins = sum(r > cpp),
      repetitions = length(r)
    )
  })
)
print(summary, row.names = FALSE)
write.csv(summary, "data-raw/bfs-native/summary.csv", row.names = FALSE)
for (name in c("parity", "validation")) {
  checks <- readRDS(paste0("data-raw/bfs-native/", name, ".rds"))
  rows <- do.call(
    rbind,
    lapply(names(checks), function(case) {
      data.frame(case = case, as.list(checks[[case]]))
    })
  )
  write.csv(
    rows,
    paste0("data-raw/bfs-native/", name, ".csv"),
    row.names = FALSE
  )
}
