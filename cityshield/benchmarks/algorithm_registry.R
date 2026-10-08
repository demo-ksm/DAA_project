# =============================================================================
# benchmarks/algorithm_registry.R - one benchmark ENTRY per algorithm
#
#   Rscript benchmarks/algorithm_registry.R
#
# Every graded algorithm is timed once at a fixed, modest size through time_it()
# and written to benchmarks/results/algorithm_registry.csv together with its
# theoretical time/space complexity (the 12 comparison experiments live in
# run_all.R; this file guarantees that EVERY algorithm has a benchmark entry).
# =============================================================================
source("benchmarks/common.R")

# input generators for the geometry entries
rand_segments_bench <- function(n, seed = 1L) {
  set.seed(seed); span <- 1000; L <- span / sqrt(n)
  x1 <- runif(n, 0, span); y1 <- runif(n, 0, span); a <- runif(n, 0, 2 * pi); len <- runif(n, 0.3 * L, L)
  segments_new(x1, y1, x1 + len * cos(a), y1 + len * sin(a))
}
hull_points <- function(n, seed = 1L) { set.seed(seed); list(x = runif(n), y = runif(n)) }

E <- function(module, algorithm, time, space, n, make, run) list(module = module, algorithm = algorithm, time = time, space = space, n = n, make = make, run = run)
rnd <- function(n, seed = 1) { set.seed(seed); runif(n) }
rtext <- function(n, alpha = letters[1:4]) { set.seed(2); utf8ToInt(paste(sample(alpha, n, TRUE), collapse = "")) }
rand_dist <- function(k) { set.seed(3); P <- cbind(runif(k), runif(k)); as.matrix(stats::dist(P)) }
und_edges <- function(n, p, seed = 5) { set.seed(seed); eu <- integer(0); ev <- integer(0)
  for (i in 1:(n - 1)) for (j in (i + 1):n) if (runif(1) < p) { eu <- c(eu, i); ev <- c(ev, j) }; list(n = n, eu = eu, ev = ev) }

registry <- list(
  # ---- Module 1
  E("M1", "fractional_knapsack", "O(n log n)", "O(n)", 2000, function(n) list(v = rnd(n), w = rnd(n, 2) + 0.1), function(x) fractional_knapsack(x$v, x$w, n_cap(x))),
  E("M1", "max_subarray_dc", "O(n log n)", "O(log n)", 2000, function(n) rnd(n) - 0.5, max_subarray_dc),
  # ---- Module 2
  E("M2", "lcs", "O(mn)", "O(mn)", 150, function(n) { set.seed(1); list(a = sample(letters[1:4], n, TRUE), b = sample(letters[1:4], n, TRUE)) }, function(x) lcs(x$a, x$b)),
  E("M2", "knapsack01_dp", "O(nW)", "O(nW)", 60, function(n) { set.seed(1); list(w = sample(1:30, n, TRUE), v = sample(1:60, n, TRUE)) }, function(x) knapsack01(x$w, x$v, 400)),
  E("M2", "tsp_bb_lc", "O(n!) worst", "O(n^2 * frontier)", 9, function(k) rand_dist(k), function(D) tsp_bb(D, "lc")),
  E("M2", "tsp_bb_lifo", "O(n!) worst", "O(n^3)", 9, function(k) rand_dist(k), function(D) tsp_bb(D, "lifo")),
  E("M2", "tsp_bb_fifo", "O(n!) worst", "O(n^2 * frontier)", 8, function(k) rand_dist(k), function(D) tsp_bb(D, "fifo")),
  # ---- Module 3
  E("M3", "kmp_match", "O(n+m)", "O(m)", 4000, function(n) rtext(n), function(t) kmp_match(t, t[100:107])),
  # ---- Module 4
  E("M4", "dijkstra", "O((V+E) log V)", "O(V+E)", 3000, function(k) sub_graph(k, 1), function(g) dijkstra(g, 1L)),
  E("M4", "push_relabel", "O(V^3)", "O(V)", 800, function(k) { g <- sub_graph(k, 3); list(fg = flow_from_csr(g, g$extra$cap), t = g$n) }, function(x) { flow_reset(x$fg); push_relabel(x$fg, 1L, x$t) }),
  E("M4", "bipartite_matching_flow", "O(min(l,r) E)", "O(V+E)", 60, function(n) { set.seed(1); k <- 4 * n; list(n = n, l = sample.int(n, k, TRUE), r = sample.int(n, k, TRUE)) }, function(x) max_bipartite_matching(x$n, x$n, x$l, x$r)),
  # ---- Module 5
  E("M5", "sweep_line", "O((n+k) log n)", "O(n+k)", 200, function(n) rand_segments_bench(n), sweep_line_intersections),
  E("M5", "graham_scan", "O(n log n)", "O(n)", 2000, function(n) hull_points(n), function(p) graham_scan(p$x, p$y)),
  # ---- Module 6
  E("M6", "randomized_quicksort", "O(n log n) exp.", "O(log n)", 2000, function(n) rnd(n), randomized_quicksort),
  E("M6", "deterministic_quicksort", "O(n^2) worst", "O(log n)", 2000, function(n) rnd(n), deterministic_quicksort),
  E("M6", "karger_run", "O(m alpha(n))", "O(n+m)", 40, function(n) { e <- csr_undirected_edges(sub_graph(n, 41)); e }, function(e) karger_run(e$n, e$eu, e$ev)),
  # ---- Module 7
  E("M7", "vertex_cover_2approx", "O(n+m)", "O(n)", 2000, function(n) und_edges(n, 3 / n), function(e) vertex_cover_approx(e$n, e$eu, e$ev)),
  E("M7", "vertex_cover_exact", "O(2^k)", "O(n)", 30, function(n) und_edges(n, 0.12), function(e) vertex_cover_exact(e$n, e$eu, e$ev)),
  E("M7", "set_cover_greedy", "O(S sum|s|)", "O(U+S)", 100, function(n) { D <- rand_dist(n); coverage_sets(D, 0.2) }, function(s) set_cover_greedy(s, length(s))),
  E("M7", "set_cover_exact", "O(f^k)", "O(U)", 40, function(n) { D <- rand_dist(n); coverage_sets(D, 0.3) }, function(s) set_cover_exact(s, length(s))),
  E("M7", "tsp_2approx", "O(n^2 log n)", "O(n)", 100, function(k) rand_dist(k), tsp_approx)
)
n_cap <- function(x) sum(x$w) / 3

rows <- NULL
for (e in registry) {
  input <- tryCatch(e$make(e$n), error = function(err) { message("make failed: ", e$algorithm, " - ", conditionMessage(err)); NULL })
  if (is.null(input)) next
  r <- tryCatch(time_it(function() e$run(input), e$algorithm, e$n, "algorithm_registry", min_iter = 2L, time_budget = 0.4),
                error = function(err) { message("run failed: ", e$algorithm, " - ", conditionMessage(err)); NULL })
  if (is.null(r)) next
  r$module <- e$module; r$theory_time <- e$time; r$theory_space <- e$space
  rows <- rbind(rows, r); cat(sprintf("%-28s n=%-6s %10.3f ms\n", e$algorithm, e$n, r$median_ms)); flush.console()
}
write_exp(rows, "algorithm_registry")
cat("\nRegistry entries:", nrow(rows), "of", length(registry), "\n")
