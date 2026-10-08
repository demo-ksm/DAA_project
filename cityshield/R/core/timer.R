# =============================================================================
# core/timer.R  -  The ONE reusable timing wrapper for every benchmark
#
# Every benchmark in the project calls time_it(); nothing else touches `bench`.
# That keeps measurements consistent (same iteration policy, same units) and
# makes the CSV files in benchmarks/results/ directly comparable.
# =============================================================================

#' Time one expression-producing function at one input size.
#'
#' @param fn        zero-argument function that runs the algorithm on a
#'                  PRE-BUILT input (build the input outside, so setup cost is
#'                  not measured).
#' @param algorithm label, e.g. "tsp_branch_bound"
#' @param n         input size (the x-axis of the benchmark)
#' @param experiment benchmark name, e.g. "maxflow_ek_vs_pr"
#' @param params    optional named list of extra columns (pattern length, ...)
#' @param min_iter,max_iter bounds on repetitions (slow algorithms run few times)
#' @param time_budget soft budget in seconds for the repetitions
#' @param memory FALSE skips allocation profiling (Rprofmem is very slow on allocation-heavy code); mem_kb is then NA
#' @param warmup FALSE skips the warm-up call (use with min_iter = 1 for runs of many seconds)
#' @return one-row data.frame: experiment, algorithm, n, median_ms, min_ms,
#'         mem_kb, iterations, plus `params`.
time_it <- function(fn, algorithm, n, experiment = "adhoc", params = list(),
                    min_iter = 3L, max_iter = 50L, time_budget = 2, warmup = TRUE, memory = TRUE) {
  if (warmup) {
    # Single warm-up call, which also tells us how slow the algorithm is.
    t0 <- proc.time()[["elapsed"]]
    fn()
    first <- proc.time()[["elapsed"]] - t0
    # Choose repetitions so the total stays near time_budget seconds.
    reps <- if (first <= 0) max_iter else floor(time_budget / first)
    reps <- max(min_iter, min(max_iter, reps))
  } else {
    reps <- min_iter                       # very slow algorithms: measure once, no warm-up
  }

  # bench::mark measures time AND allocated memory; check = FALSE because we
  # compare algorithms elsewhere (in tests), not here.
  res <- suppressWarnings(
    bench::mark(fn(), iterations = reps, check = FALSE, memory = memory,
                time_unit = "ms", filter_gc = FALSE))
  row <- data.frame(
    experiment = experiment,
    algorithm  = algorithm,
    n          = n,
    median_ms  = as.numeric(res$median),
    min_ms     = as.numeric(res$min),
    mem_kb     = if (memory) as.numeric(res$mem_alloc) / 1024 else NA_real_,
    iterations = as.integer(res$n_itr),
    stringsAsFactors = FALSE
  )
  for (nm in names(params)) row[[nm]] <- params[[nm]]
  row
}

#' Append benchmark rows to benchmarks/results/<experiment>.csv.
#' Creates the file (with header) if needed. `results_dir` defaults to the
#' project's benchmarks/results folder.
save_bench <- function(rows, experiment, results_dir = file.path(PROJECT_ROOT, "benchmarks", "results"),
                       overwrite = FALSE) {
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  path <- file.path(results_dir, paste0(experiment, ".csv"))
  if (overwrite || !file.exists(path)) {
    utils::write.csv(rows, path, row.names = FALSE)
  } else {
    old <- utils::read.csv(path, stringsAsFactors = FALSE)
    all_cols <- union(names(old), names(rows))
    for (cn in setdiff(all_cols, names(old))) old[[cn]] <- NA
    for (cn in setdiff(all_cols, names(rows))) rows[[cn]] <- NA
    utils::write.csv(rbind(old[, all_cols], rows[, all_cols]), path, row.names = FALSE)
  }
  invisible(path)
}

#' Run time_it() over a vector of sizes and bind the rows together.
#'
#' @param make_input function(n) -> input object (built OUTSIDE the timed part)
#' @param run        function(input) -> runs the algorithm
#' @param sizes      numeric vector of input sizes
#' @param stop_above_s stop growing sizes once a single run exceeds this many
#'                   seconds (documents the "practical limit" in R)
bench_sizes <- function(make_input, run, sizes, algorithm, experiment,
                        params = list(), stop_above_s = 10, ...) {
  out <- NULL
  for (n in sizes) {
    input <- make_input(n)
    row <- time_it(function() run(input), algorithm, n, experiment, params, ...)
    out <- if (is.null(out)) row else rbind(out, row)
    if (row$median_ms / 1000 > stop_above_s) {
      message(sprintf("[%s] %s stopped at n=%s (%.1fs per run > cap)",
                      experiment, algorithm, n, row$median_ms / 1000))
      break
    }
  }
  out
}
