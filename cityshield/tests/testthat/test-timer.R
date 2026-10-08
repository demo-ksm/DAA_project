test_that("time_it returns one well-formed row", {
  r <- time_it(function() { s <- 0; for (i in 1:1000) s <- s + i; s },
               algorithm = "sum_loop", n = 1000, experiment = "unit_test",
               params = list(note = "x"), time_budget = 0.2)
  expect_equal(nrow(r), 1L)
  expect_true(all(c("experiment", "algorithm", "n", "median_ms", "min_ms",
                    "mem_kb", "iterations", "note") %in% names(r)))
  expect_gt(r$median_ms, 0)
  expect_gte(r$iterations, 3L)
})

test_that("a larger input is measured as slower", {
  work <- function(n) function() { s <- 0; for (i in seq_len(n)) s <- s + i; s }
  small <- time_it(work(1e3), "loop", 1e3, time_budget = 0.2)
  large <- time_it(work(1e5), "loop", 1e5, time_budget = 0.2)
  expect_gt(large$median_ms, small$median_ms)
})

test_that("save_bench creates then appends to a CSV", {
  d <- tempfile(); dir.create(d)
  r1 <- data.frame(experiment = "e", algorithm = "a", n = 1, median_ms = 1)
  r2 <- data.frame(experiment = "e", algorithm = "b", n = 2, median_ms = 2, extra = "z")
  p <- save_bench(r1, "e", results_dir = d)
  save_bench(r2, "e", results_dir = d)
  back <- read.csv(p)
  expect_equal(nrow(back), 2L)
  expect_true("extra" %in% names(back))
  save_bench(r1, "e", results_dir = d, overwrite = TRUE)
  expect_equal(nrow(read.csv(p)), 1L)
})

test_that("bench_sizes stops once a run exceeds the cap", {
  slow <- function(n) { t0 <- proc.time()[["elapsed"]]; while (proc.time()[["elapsed"]] - t0 < n / 1000) NULL }
  res <- suppressMessages(bench_sizes(function(n) n, slow, sizes = c(10, 20, 1500, 2000),
                                      algorithm = "busy", experiment = "t",
                                      stop_above_s = 1, min_iter = 1L, max_iter = 1L, time_budget = 0))
  expect_equal(res$n, c(10, 20, 1500))
})
