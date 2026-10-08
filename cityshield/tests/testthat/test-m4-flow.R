# Validation of push-relabel max flow, bipartite matching and the evacuation glue.

flow_value_pr <- function(fe, s, t, ...) push_relabel(build_fg(fe), s, t, ...)$value

test_that("CLRS textbook network has max flow 23", {
  fe <- list(n = 6, from = c(1, 1, 2, 3, 2, 4, 3, 5, 4, 5),
             to = c(2, 3, 4, 2, 3, 3, 5, 4, 6, 6),
             cap = c(16, 13, 12, 4, 10, 9, 14, 7, 20, 4))
  fe$from <- as.integer(fe$from); fe$to <- as.integer(fe$to)
  expect_equal(flow_value_pr(fe, 1L, 6L), 23)
})

test_that("push-relabel agrees with igraph::max_flow on random networks", {
  for (seed in 1:8) {
    fe <- rand_flow_edges(25, 90, seed)
    ig <- igraph::make_graph(rbind(fe$from, fe$to), n = fe$n, directed = TRUE)
    ref <- igraph::max_flow(ig, 1, 25, capacity = fe$cap)$value
    expect_equal(flow_value_pr(fe, 1L, 25L), ref, info = paste("seed", seed))
  }
})

test_that("flow is valid: capacity respected, conservation holds, max-flow = min-cut", {
  fe <- rand_flow_edges(20, 80, 3)
  fg <- build_fg(fe)
  val <- push_relabel(fg, 1L, 20L, complete = TRUE)$value
  f <- flow_values(fg)
  net <- numeric(fe$n)
  for (k in seq_along(fe$from)) {
    e <- 2L * k - 1L
    expect_gte(f[e], -1e-9); expect_lte(f[e], fe$cap[k] + 1e-9)
    net[fe$from[k]] <- net[fe$from[k]] - f[e]
    net[fe$to[k]]   <- net[fe$to[k]] + f[e]
  }
  for (v in 2:19) expect_equal(net[v], 0, info = paste("node", v))   # conservation
  expect_equal(net[20], val)
  side <- min_cut_side(fg, 1L)                       # min cut capacity == flow value
  cutcap <- 0
  for (k in seq_along(fe$from)) if (side[fe$from[k]] && !side[fe$to[k]]) cutcap <- cutcap + fe$cap[k]
  expect_equal(cutcap, val)
})

test_that("flow edge cases: no path, single edge, parallel edges, zero capacity, reset", {
  none <- list(n = 3, from = 1L, to = 2L, cap = 5L)
  expect_equal(flow_value_pr(none, 1L, 3L), 0)
  one <- list(n = 2, from = 1L, to = 2L, cap = 7L)
  expect_equal(flow_value_pr(one, 1L, 2L), 7)
  par <- list(n = 2, from = c(1L, 1L, 1L), to = c(2L, 2L, 2L), cap = c(1L, 2L, 3L))
  expect_equal(flow_value_pr(par, 1L, 2L), 6)
  zero <- list(n = 3, from = c(1L, 2L), to = c(2L, 3L), cap = c(0L, 4L))
  expect_equal(flow_value_pr(zero, 1L, 3L), 0)
  bott <- list(n = 4, from = c(1L, 2L, 3L), to = c(2L, 3L, 4L), cap = c(10L, 1L, 10L))
  expect_equal(flow_value_pr(bott, 1L, 4L), 1)
  fg <- build_fg(one); push_relabel(fg, 1L, 2L); flow_reset(fg)
  expect_equal(push_relabel(fg, 1L, 2L)$value, 7)               # a second run after reset gives the same answer
})

test_that("push-relabel: phase-1 value = complete flow value, with or without gap / global relabel", {
  for (seed in 1:8) {
    fe <- rand_flow_edges(25, 70, seed)
    ig <- igraph::make_graph(rbind(fe$from, fe$to), n = fe$n, directed = TRUE)
    ref <- igraph::max_flow(ig, 1, 25, capacity = fe$cap)$value
    for (gr in c(TRUE, FALSE)) for (gp in c(TRUE, FALSE)) for (cm in c(TRUE, FALSE))
      expect_equal(push_relabel(build_fg(fe), 1L, 25L, global_relabel = gr, gap = gp, complete = cm)$value, ref,
                   info = paste(seed, gr, gp, cm))
  }
})

test_that("push-relabel phase 1 is fast when most excess cannot reach the sink", {
  # a wide fan of dead-end nodes hanging off the source: phase 1 must not shuffle their excess around
  n <- 300; fg <- flow_new(n, 2 * n)
  for (v in 2:(n - 2)) { flow_add_edge(fg, 1L, v, 10); flow_add_edge(fg, v, (v %% (n - 3)) + 2L, 10) }
  flow_add_edge(fg, 2L, n, 5)
  r <- push_relabel(fg, 1L, n)
  expect_equal(r$value, 5)
  fg2 <- flow_reset(fg); r2 <- push_relabel(fg2, 1L, n, complete = TRUE)
  expect_equal(r2$value, 5)
  expect_lt(r$relabels, r2$relabels + 1)
})

# ---- bipartite matching -----------------------------------------------------

test_that("matching agrees with igraph and brute force", {
  for (seed in 1:10) {
    set.seed(seed)
    nl <- sample(3:7, 1); nr <- sample(3:7, 1); m <- sample(4:18, 1)
    pairs <- unique(cbind(sample.int(nl, m, TRUE), sample.int(nr, m, TRUE)))
    l <- pairs[, 1]; r <- pairs[, 2]
    ig <- igraph::make_bipartite_graph(c(rep(FALSE, nl), rep(TRUE, nr)),
                                       as.vector(rbind(l, nl + r)))
    ref <- igraph::max_bipartite_match(ig)$matching_size
    {
      res <- max_bipartite_matching(nl, nr, l, r)
      expect_equal(res$size, ref, info = paste("push-relabel", seed))
      # validity: no shared endpoints, only allowed pairs
      expect_equal(anyDuplicated(res$pairs[, 1]), 0L); expect_equal(anyDuplicated(res$pairs[, 2]), 0L)
      for (i in seq_len(nrow(res$pairs))) {
        ok <- FALSE
        for (k in seq_along(l)) if (l[k] == res$pairs[i, 1] && r[k] == res$pairs[i, 2]) ok <- TRUE
        expect_true(ok)
      }
    }
    expect_equal(matching_brute(nl, nr, l, r), ref, info = paste("brute", seed))
  }
})

test_that("matching edge cases", {
  expect_equal(max_bipartite_matching(3, 3, integer(0), integer(0))$size, 0L)
  expect_equal(max_bipartite_matching(1, 1, 1L, 1L)$size, 1L)
  # complete 3x3 -> perfect matching
  l <- rep(1:3, each = 3); r <- rep(1:3, 3)
  expect_equal(max_bipartite_matching(3, 3, l, r)$size, 3L)
  # star: one right node wanted by everyone -> only one match
  expect_equal(max_bipartite_matching(4, 1, 1:4, rep(1L, 4))$size, 1L)
  expect_equal(matching_brute(4, 1, 1:4, rep(1L, 4)), 1L)
  # augmenting-path needed: greedy would match 1-1 and block 2
  expect_equal(max_bipartite_matching(2, 2, c(1L, 1L, 2L), c(1L, 2L, 1L))$size, 2L)
})

test_that("evacuation network moves people as far as roads/shelters allow", {
  # line 1 -> 2 -> 3 with capacity 2 veh/h, 3 persons per vehicle => 6 per edge
  g <- csr_build(3, c(1L, 2L), c(2L, 3L), c(10, 10), extra = list(cap = c(2, 2)))
  net <- evacuation_network(g, danger = 1L, people = 100, shelters = 3L, shelter_cap = 50)
  expect_equal(push_relabel(net$fg, net$s, net$t, complete = TRUE)$value, 6)       # road bottleneck
  ef <- evacuation_edge_flow(net)
  expect_equal(ef, c(6, 6))
  net2 <- evacuation_network(g, 1L, 4, 3L, 50)
  expect_equal(push_relabel(net2$fg, net2$s, net2$t)$value, 4)    # supply bottleneck
})

test_that("dispatch_feasibility respects the time limit", {
  g <- csr_build(4, c(1L, 2L, 3L), c(2L, 3L, 4L), c(5, 5, 5))
  f <- dispatch_feasibility(g, amb_nodes = 1L, inc_nodes = c(2L, 3L, 4L), max_time = 10)
  expect_equal(f$r, c(1L, 2L)); expect_equal(f$time, c(5, 10))
})
