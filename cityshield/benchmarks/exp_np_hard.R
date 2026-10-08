# =============================================================================
# Benchmarks 5-8: the NP-hard exact-vs-approximate study + Karger
#   5 TSP branch & bound vs 2-approximation   6 exact vs greedy set cover
#   7 exact vs 2-approx vertex cover          8 Karger success rate vs trials
# Instances are 50-150 node road subgraphs (or sites/universes drawn from them).
# =============================================================================

# ---- 5. TSP: B&B vs 2-approximation -----------------------------------------------------------
exp05_tsp <- function(quick = FALSE) {
  ks <- if (quick) 5:7 else 5:13
  reps <- if (quick) 2 else 4
  node_cap <- if (quick) 3000 else 40000
  base <- connected_instance(if (quick) 40 else 100, 11L)
  out <- NULL
  for (k in ks) for (r in seq_len(reps)) {
    set.seed(1000 * k + r)
    sites <- sample.int(base$g$n, k)
    D <- base$D[sites, sites]
    ap <- tsp_approx(D)
    t_ap <- time_it(function() tsp_approx(D), "tsp_2approx", k, "tsp", min_iter = 3L, time_budget = 0.3)
    res <- NULL
    t_bb <- time_it(function() { res <<- tsp_bb(D, "lc", max_nodes = node_cap) }, "tsp_branch_bound", k, "tsp",
                    min_iter = 1L, warmup = FALSE, memory = FALSE)
    opt <- if (res$complete) res$cost else NA_real_          # ratio only where B&B proved the optimum
    out <- rbind(out, data.frame(n = k, rep = r, bb_ms = t_bb$median_ms, approx_ms = t_ap$median_ms,
                                 bb_complete = res$complete, bb_expanded = res$expanded, bb_generated = res$generated,
                                 optimal = opt, approx_cost = ap$cost, mst_weight = ap$mst_weight,
                                 ratio = ap$cost / opt, graph = bench_graph_source()))
  }
  write_exp(out, "05_tsp_bb_vs_approx")
  agg <- do.call(rbind, lapply(split(out, out$n), function(d)
    data.frame(n = d$n[1], bb = med(d$bb_ms), approx = med(d$approx_ms), expanded = med(d$bb_expanded))))
  long <- rbind(data.frame(n = agg$n, algorithm = "branch_and_bound", median_ms = agg$bb),
                data.frame(n = agg$n, algorithm = "2_approximation", median_ms = agg$approx))
  p <- runtime_plot(long, "TSP inspection routes: exact B&B vs 2-approximation",
                    "Median over random inspection-site sets; B&B grows super-polynomially, the approximation stays flat.",
                    xlab = "number of inspection sites k")
  save_plot(p, "05_tsp_runtime")
  p2 <- ggplot(out, aes(factor(n), ratio)) +
    geom_hline(yintercept = 2, colour = "#d03b3b", linetype = "dashed") +
    geom_boxplot(fill = "#cde2fb", colour = "#2a78d6", outlier.colour = "#52514e") +
    labs(title = "Approximation ratio of the MST-based tour", subtitle = "dashed line = theoretical bound 2",
         x = "number of inspection sites k", y = "tour cost / optimal cost") + bench_theme()
  save_plot(p2, "05_tsp_ratio")
  limit <- agg$n[agg$bb > 5000]
  list(data = out)
}

# ---- 6. set cover: exact vs greedy ---------------------------------------------------------------------
exp06_setcover <- function(quick = FALSE) {
  sizes <- if (quick) c(15, 25, 35) else c(20, 30, 40, 50, 60, 80, 100)
  reps <- if (quick) 2 else 5
  node_cap <- if (quick) 3000 else 60000
  base <- connected_instance(if (quick) 40 else 110, 21L)
  out <- NULL
  for (U in sizes) for (r in seq_len(reps)) {
    set.seed(2000 * U + r)
    idx <- sample.int(base$g$n, U)
    D <- base$D[idx, idx]
    rad <- stats::quantile(D[upper.tri(D)], 0.12 + 0.03 * (r - 1))     # varying coverage radius
    sets <- coverage_sets(D, rad)
    gr <- set_cover_greedy(sets, U)
    ex <- NULL
    t_ex <- time_it(function() { ex <<- set_cover_exact(sets, U, max_nodes = node_cap) }, "set_cover_exact", U, "setcover",
                    min_iter = 1L, warmup = FALSE, memory = FALSE)
    t_gr <- time_it(function() set_cover_greedy(sets, U), "set_cover_greedy", U, "setcover", min_iter = 3L, time_budget = 0.3)
    d <- max(lengths(sets))
    out <- rbind(out, data.frame(n = U, rep = r, radius = as.numeric(rad), max_set = d, exact_ms = t_ex$median_ms,
                                 greedy_ms = t_gr$median_ms, exact_size = ex$size, greedy_size = gr$size,
                                 exact_complete = ex$complete, ratio = gr$size / ex$size,
                                 harmonic_bound = harmonic_bound(d), ln_bound = log(U) + 1, graph = bench_graph_source()))
  }
  write_exp(out, "06_set_cover")
  agg <- do.call(rbind, lapply(split(out, out$n), function(d)
    data.frame(n = d$n[1], exact = med(d$exact_ms), greedy = med(d$greedy_ms))))
  long <- rbind(data.frame(n = agg$n, algorithm = "exact_branch_bound", median_ms = agg$exact),
                data.frame(n = agg$n, algorithm = "greedy", median_ms = agg$greedy))
  p <- runtime_plot(long, "Station placement (set cover): exact vs greedy",
                    "Universe = demand nodes of a road subgraph; stations cover nodes within a travel-time radius.",
                    xlab = "universe size n")
  save_plot(p, "06_set_cover_runtime")
  p2 <- ggplot(out, aes(n, ratio)) +
    geom_line(aes(y = harmonic_bound, group = rep), colour = "#c3c2b7") +
    geom_point(colour = "#2a78d6", size = 2) +
    scale_y_continuous(limits = c(1, NA)) +
    labs(title = "Greedy set cover quality", subtitle = "points = greedy/optimal; grey = H(d) bound (d = largest set)",
         x = "universe size n", y = "greedy size / optimal size") + bench_theme()
  save_plot(p2, "06_set_cover_ratio")
  ok <- out[out$exact_complete, ]
  list(data = out)
}

# ---- 7. vertex cover: exact vs 2-approximation ---------------------------------------------------------------
exp07_vertexcover <- function(quick = FALSE) {
  sizes <- if (quick) c(10, 16, 22) else c(10, 20, 30, 40, 50, 60, 80)
  node_cap <- if (quick) 5000 else 150000
  out <- NULL
  for (n in sizes) for (fam in c("road", "random", "matching")) {
    for (r in 1:(if (quick) 2 else 4)) {
      set.seed(3000 * n + r)
      if (fam == "road") {
        e <- csr_undirected_edges(sub_graph(n, 30L + r)); nn <- e$n; eu <- e$eu; ev <- e$ev
      } else if (fam == "random") {
        nn <- n; eu <- integer(0); ev <- integer(0)
        for (i in 1:(n - 1)) for (j in (i + 1):n) if (runif(1) < 2.5 / n) { eu <- c(eu, i); ev <- c(ev, j) }
      } else {                                              # perfect matching: approx is exactly 2x optimal
        nn <- n - n %% 2; eu <- seq(1L, nn - 1L, by = 2L); ev <- eu + 1L
      }
      if (length(eu) == 0L) next
      ex <- NULL
      t_ex <- time_it(function() { ex <<- vertex_cover_exact(nn, eu, ev, max_nodes = node_cap) }, "vertex_cover_exact", n, "vc",
                      min_iter = 1L, warmup = FALSE, memory = FALSE)
      t_ap <- time_it(function() vertex_cover_approx(nn, eu, ev), "vertex_cover_2approx", n, "vc", min_iter = 3L, time_budget = 0.2)
      ap <- vertex_cover_approx(nn, eu, ev)
      out <- rbind(out, data.frame(n = n, family = fam, rep = r, nodes = nn, edges = length(eu),
                                   exact_ms = t_ex$median_ms, approx_ms = t_ap$median_ms, exact_size = ex$size,
                                   approx_size = ap$size, exact_complete = ex$complete, ratio = ap$size / ex$size,
                                   bnb_nodes = ex$nodes, graph = if (fam == "road") bench_graph_source() else "synthetic"))
    }
  }
  write_exp(out, "07_vertex_cover")
  agg <- do.call(rbind, lapply(split(out, list(out$n, out$family)), function(d)
    data.frame(n = d$n[1], family = d$family[1], exact = med(d$exact_ms), approx = med(d$approx_ms))))
  long <- rbind(data.frame(n = agg$n, family = agg$family, algorithm = "exact_branch_bound", median_ms = agg$exact),
                data.frame(n = agg$n, family = agg$family, algorithm = "2_approximation", median_ms = agg$approx))
  p <- runtime_plot(long, "Sensor placement (vertex cover): exact vs 2-approximation",
                    "Road subgraphs, sparse random graphs, and perfect matchings (the approximation's worst case).",
                    facet = "family", xlab = "number of nodes")
  save_plot(p, "07_vertex_cover_runtime", 9, 4.2)
  p2 <- ggplot(out[out$exact_complete, ], aes(n, ratio, colour = family)) +
    geom_hline(yintercept = 2, colour = "#d03b3b", linetype = "dashed") +
    geom_point(size = 2.2, alpha = 0.85) + scale_colour_manual(values = SERIES) +
    scale_y_continuous(limits = c(1, 2.15)) +
    labs(title = "2-approximation quality for vertex cover", subtitle = "dashed = guaranteed bound 2",
         x = "number of nodes", y = "approx size / optimal size") + bench_theme()
  save_plot(p2, "07_vertex_cover_ratio")
  ok <- out[out$exact_complete, ]
  list(data = out)
}

# ---- 8. Karger: success rate vs number of trials ---------------------------------------------------------------
karger_instance <- function(family, n, seed) {
  if (family == "road") {
    e <- csr_undirected_edges(sub_graph(n, seed)); list(n = e$n, eu = e$eu, ev = e$ev)
  } else if (family == "two_cliques") {                    # two K_{n/2} joined by 2 edges: min cut 2, many near-cuts
    a <- n %/% 2; eu <- integer(0); ev <- integer(0)
    for (blk in 0:1) for (i in 1:(a - 1)) for (j in (i + 1):a) { eu <- c(eu, blk * a + i); ev <- c(ev, blk * a + j) }
    list(n = 2 * a, eu = c(eu, 1L, 2L), ev = c(ev, a + 1L, a + 2L))
  } else {                                                 # cycle: n(n-1)/2 distinct minimum cuts of size 2 (easy)
    list(n = n, eu = 1:n, ev = c(2:n, 1L))
  }
}

exp08_karger <- function(quick = FALSE) {
  fams <- c("road", "two_cliques", "cycle")
  sizes <- if (quick) c(12, 20) else c(20, 40, 80)
  runs <- if (quick) 150 else 600
  Ts <- c(1, 2, 5, 10, 20, 50, 100, 200, 400)
  out <- NULL; meta <- NULL
  for (fam in fams) for (n in sizes) {
    inst <- karger_instance(fam, n, 41L)
    if (!edges_connected(inst$n, inst$eu, inst$ev)) next
    ig <- igraph::make_graph(as.vector(rbind(inst$eu, inst$ev)), n = inst$n, directed = FALSE)
    true_cut <- igraph::min_cut(ig)                                   # validation (allowed: igraph)
    set.seed(7)
    hit <- logical(runs); cuts <- numeric(runs)
    t_run <- time_it(function() karger_run(inst$n, inst$eu, inst$ev), "karger_run", inst$n, "karger",
                     min_iter = 3L, time_budget = 0.5)
    for (i in seq_len(runs)) { cuts[i] <- karger_run(inst$n, inst$eu, inst$ev)$cut; hit[i] <- cuts[i] == true_cut }
    p1 <- mean(hit); bound <- 2 / (inst$n * (inst$n - 1))
    for (Tt in Ts) {                                                  # success after T trials: resample groups of T runs
      reps <- 400; ok <- 0
      for (r in seq_len(reps)) if (any(hit[sample.int(runs, Tt, replace = TRUE)])) ok <- ok + 1
      out <- rbind(out, data.frame(family = fam, n = inst$n, trials = Tt, empirical = ok / reps,
                                   from_p1 = 1 - (1 - p1)^Tt, theory_bound = 1 - (1 - bound)^Tt))
    }
    meta <- rbind(meta, data.frame(family = fam, n = inst$n, edges = length(inst$eu), min_cut_igraph = true_cut,
                                   karger_min_found = min(cuts), per_trial_success = p1, guaranteed_bound = bound,
                                   run_ms = t_run$median_ms, trials_for_99pct = karger_trials_needed(inst$n, 0.01),
                                   graph = if (fam == "road") bench_graph_source() else "synthetic"))
  }
  write_exp(out, "08_karger_success")
  write_exp(meta, "08_karger_meta")
  out$case <- paste0(out$family, " (n=", out$n, ")")
  p <- ggplot(out, aes(trials)) +
    geom_line(aes(y = theory_bound, colour = "theoretical bound 1-(1-2/n(n-1))^T"), linewidth = 0.7, linetype = "dashed") +
    geom_line(aes(y = empirical, colour = "observed success rate"), linewidth = 0.9) +
    geom_point(aes(y = empirical, colour = "observed success rate"), size = 1.6) +
    scale_x_log10() + scale_colour_manual(values = c("observed success rate" = "#2a78d6",
                                                     "theoretical bound 1-(1-2/n(n-1))^T" = "#898781")) +
    facet_wrap(~ case, ncol = 3) +
    labs(title = "Karger's min cut: success probability vs number of trials",
         subtitle = "observed always sits above the guaranteed bound; the bound is only tight for cut-rich dense graphs",
         x = "independent trials T (log scale)", y = "P(min cut found)") + bench_theme()
  save_plot(p, "08_karger_success", 10, 7)
  list(data = out)
}
