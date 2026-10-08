# =============================================================================
# Benchmark 9: quicksort (randomized vs deterministic pivot on adversarial inputs).
# =============================================================================

# ---- 9. quicksort ------------------------------------------------------------------------
exp09_quicksort <- function(quick = FALSE) {
  sizes <- if (quick) c(100, 200, 400) else c(250, 500, 1000, 2000, 4000)
  kinds <- c("random", "sorted", "reversed", "organ")
  out <- NULL
  for (kind in kinds) for (alg in c("randomized", "deterministic")) {
    f <- if (alg == "randomized") randomized_quicksort else deterministic_quicksort
    r <- bench_sizes(function(n) { set.seed(5); quicksort_input(n, kind) }, function(x) f(x), sizes, alg, "quicksort",
                     params = list(input = kind), stop_above_s = 8, min_iter = 1L, time_budget = 1)
    for (i in seq_len(nrow(r))) {                      # comparisons (mean of 5 runs for the randomized pivot)
      set.seed(5); x <- quicksort_input(r$n[i], kind)
      r$comparisons[i] <- mean(vapply_num(5L, function(j) f(x)$comparisons))
    }
    out <- rbind(out, r)
  }
  write_exp(out, "09_quicksort")
  out$label <- paste(out$algorithm, out$input)
  p <- ggplot(out, aes(n, comparisons, colour = algorithm, group = algorithm)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) + facet_wrap(~ input, scales = "free_y") +
    scale_colour_manual(values = SERIES) + scale_x_log10() + scale_y_log10() +
    labs(title = "Quicksort comparisons: randomized vs deterministic (last-element pivot)",
         subtitle = "sorted / reversed / organ-pipe inputs are adversarial for the fixed pivot (n^2); random pivot stays ~ 2 n ln n",
         x = "input size n", y = "comparisons (log scale)") + bench_theme()
  save_plot(p, "09_quicksort_comparisons", 9, 5.5)
  save_plot(runtime_plot(out, "Quicksort runtime", facet = "input"), "09_quicksort_runtime", 9, 5.5)
  sl <- function(a, k, col = "comparisons") { d <- out[out$algorithm == a & out$input == k, ]; loglog_slope(d$n, d[[col]]) }
  list(data = out)
}

# small helper: apply f to 1..k and return a numeric vector (explicit loop; benchmark code only)
vapply_num <- function(k, f) { o <- numeric(k); for (j in seq_len(k)) o[j] <- f(j); o }
