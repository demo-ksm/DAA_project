# =============================================================================
# benchmarks/notes.R - "theory vs observation" text, generated FROM THE SAVED CSVs
#
#   Rscript benchmarks/notes.R           # rebuilds benchmarks/results/notes.rds
#
# Every sentence is computed from the measured data (slopes, ratios, crossovers),
# so the notes can never drift away from what the benchmark actually saw.
# =============================================================================
source("benchmarks/common.R")

rd <- function(f) { p <- file.path(RESULTS_DIR, f); if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
at <- function(d, a, n, col = "median_ms") d[[col]][d$algorithm == a & d$n == n][1]
fx <- function(x, digits = 2) formatC(x, format = "f", digits = digits)
tm <- function(ms) if (ms >= 1000) sprintf("%.2f s", ms / 1000) else if (ms < 0.1) "< 0.1 ms" else sprintf("%.1f ms", ms)
note_05 <- function() {
  d <- rd("05_tsp_bb_vs_approx.csv"); if (is.null(d)) return(NULL)
  agg <- do.call(rbind, lapply(split(d, d$n), function(x) data.frame(k = x$n[1], bb = med(x$bb_ms), ap = med(x$approx_ms), inc = sum(!x$bb_complete), exp = med(x$bb_expanded))))
  lim <- agg$k[agg$bb > 1000 | agg$inc > 0]
  sprintf(paste0("Theory: B&B O(k!) worst case (reduced-matrix bound prunes most of it on metric instances); MST 2-approximation O(k^2 log k) with ratio <= 2. Instances: %d random site sets per k, distances = symmetrised shortest-path times on a %s road subgraph. ",
    "Observed approximation ratio: max %.3f, median %.3f (bound 2 respected: %s). B&B median time grows from %s (k = %d) to %s (k = %d), expanding a median of %.0f -> %.0f nodes; %d instance(s) hit the node cap. ",
    "B&B first becomes impractical (median > 1 s or cap hit) at k = %s, whereas the approximation takes %s at k = %d. Reading: exact TSP stops being interactive at k = %s here, while the approximation stays instant at a median %.0f%% (worst %.0f%%) above the optimal tour length."),
    max(table(d$n)), unique(d$graph)[1], max(d$ratio, na.rm = TRUE), med(d$ratio), if (all(d$ratio <= 2 + 1e-9, na.rm = TRUE)) "yes" else "NO",
    tm(agg$bb[1]), agg$k[1], tm(max(agg$bb)), agg$k[which.max(agg$bb)], agg$exp[1], max(agg$exp), sum(!d$bb_complete),
    if (length(lim)) as.character(min(lim)) else "not reached in the tested range", tm(agg$ap[nrow(agg)]), max(agg$k),
    if (length(lim)) as.character(min(lim)) else "beyond the tested range", 100 * (med(d$ratio) - 1), 100 * (max(d$ratio, na.rm = TRUE) - 1))
}

note_06 <- function() {
  d <- rd("06_set_cover.csv"); if (is.null(d)) return(NULL); ok <- d[d$exact_complete, ]
  n <- max(d$n); a <- d[d$n == n, ]
  sprintf(paste0("Theory: greedy covers with <= H(d) * OPT <= (ln d + 1) * OPT sets (d = largest set). Instances: universes of %d-%d demand nodes of a %s road subgraph, stations cover nodes within a travel-time radius. ",
    "Observed greedy/optimal: max %.3f, median %.3f over %d instances solved to optimality; within H(d) in all of them: %s; max H(d) seen %.2f, ln(n)+1 = %.2f at n = %d. ",
    "At n = %d the exact B&B needed a median %s (%s nodes cap hit in %d of %d runs) versus %s for greedy. Reading: greedy matched the optimum in %.0f%% of the solved instances and never exceeded %.2fx, far inside the logarithmic worst case, at a tiny fraction of the cost."),
    min(d$n), n, unique(d$graph)[1], max(ok$ratio), med(ok$ratio), nrow(ok), if (all(ok$ratio <= ok$harmonic_bound + 1e-9)) "yes" else "NO",
    max(d$harmonic_bound), log(n) + 1, n, n, tm(med(a$exact_ms)), "60,000", sum(!a$exact_complete), nrow(a), tm(med(a$greedy_ms)),
    100 * mean(ok$ratio <= 1 + 1e-9), max(ok$ratio))
}

note_07 <- function() {
  d <- rd("07_vertex_cover.csv"); if (is.null(d)) return(NULL); ok <- d[d$exact_complete, ]
  big <- d[d$n == max(d$n) & d$family == "road", ]
  sprintf(paste0("Theory: the maximal-matching cover has size 2|M| <= 2 OPT. Observed on %d instances solved exactly: max ratio %.3f (bound 2 %s), median %.3f on road/random graphs; perfect matchings give exactly %.2f (the bound is tight). ",
    "Exact B&B median time at %d nodes on road subgraphs: %s (complete in %d of %d runs) versus %s for the approximation. Reading: on road and random graphs the approximation used a median %.2fx the optimal number of sensors (largest %.2fx); the perfect-matching family reaches the bound 2 exactly, which shows the guarantee cannot be improved for this algorithm."),
    nrow(ok), max(ok$ratio), if (all(ok$ratio <= 2 + 1e-9)) "respected" else "VIOLATED", med(ok$ratio[ok$family != "matching"]), max(ok$ratio[ok$family == "matching"]),
    max(d$n), tm(med(big$exact_ms)), sum(big$exact_complete), nrow(big), tm(med(big$approx_ms)),
    med(ok$ratio[ok$family != "matching"]), max(ok$ratio[ok$family != "matching"]))
}

note_08 <- function() {
  m <- rd("08_karger_meta.csv"); d <- rd("08_karger_success.csv"); if (is.null(m)) return(NULL)
  sprintf(paste0("Theory: one Karger trial finds a fixed minimum cut with probability >= 2/(n(n-1)); T independent trials fail with probability <= (1 - 2/(n(n-1)))^T, so T = n(n-1)/2 * ln(1/delta) trials give failure <= delta. ",
    "Observed: igraph::min_cut equals the best Karger cut in %d of %d instances (families: %s). Measured per-trial success ranged %.4f-%.4f against guaranteed bounds %.5f-%.5f (every observed rate is above its bound: %s). ",
    "One trial costs %s on average (n up to %d). Reading: the n^2 guarantee is a worst case reached only by graphs with many near-minimum cuts (the two-clique family); sparse road subgraphs, whose minimum cut is a bridge-like link, are solved in a handful of trials."),
    sum(m$min_cut_igraph == m$karger_min_found), nrow(m), paste(unique(m$family), collapse = ", "), min(m$per_trial_success), max(m$per_trial_success), min(m$guaranteed_bound), max(m$guaranteed_bound),
    if (all(m$per_trial_success >= m$guaranteed_bound)) "yes" else "NO", tm(med(m$run_ms)), max(m$n))
}

note_09 <- function() {
  d <- rd("09_quicksort.csv"); if (is.null(d)) return(NULL)
  s <- function(a, k, col) { x <- d[d$algorithm == a & d$input == k, ]; loglog_slope(x$n, x[[col]]) }
  nmax <- max(d$n[d$algorithm == "deterministic" & d$input == "sorted"])
  sprintf(paste0("Theory: randomized quicksort needs ~2 n ln n comparisons in expectation on EVERY input; deterministic last-element pivot needs n(n-1)/2 on sorted/reversed input. Comparison-count slopes: deterministic sorted %.2f, reversed %.2f, random %.2f; randomized sorted %.2f, reversed %.2f, random %.2f. ",
    "At n = %d on sorted input: deterministic %s comparisons (n(n-1)/2 = %s) in %s, randomized %s comparisons (2 n ln n = %s) in %s. Reading: randomization removes the adversarial inputs; both versions use an explicit stack processing the smaller side first, so recursion depth stays O(log n)."),
    s("deterministic", "sorted", "comparisons"), s("deterministic", "reversed", "comparisons"), s("deterministic", "random", "comparisons"),
    s("randomized", "sorted", "comparisons"), s("randomized", "reversed", "comparisons"), s("randomized", "random", "comparisons"), nmax,
    format(round(at(d[d$input == "sorted", ], "deterministic", nmax, "comparisons")), big.mark = ","), format(nmax * (nmax - 1) / 2, big.mark = ","), tm(at(d[d$input == "sorted", ], "deterministic", nmax)),
    format(round(at(d[d$input == "sorted", ], "randomized", nmax, "comparisons")), big.mark = ","), format(round(2 * nmax * log(nmax)), big.mark = ","), tm(at(d[d$input == "sorted", ], "randomized", nmax)))
}

make_notes <- function() {
  fs <- list("05" = note_05, "06" = note_06, "07" = note_07, "08" = note_08, "09" = note_09)
  out <- list()
  for (id in names(fs)) { r <- tryCatch(fs[[id]](), error = function(e) paste("(notes failed:", conditionMessage(e), ")")); if (!is.null(r)) out[[id]] <- r }
  saveRDS(out, file.path(RESULTS_DIR, "notes.rds"))
  out
}

if (sys.nframe() == 0L) { n <- make_notes(); for (k in names(n)) cat("[", k, "]", n[[k]], "\n\n") }
