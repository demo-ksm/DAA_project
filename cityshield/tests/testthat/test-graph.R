# Random directed test graph (as plain edge lists).
make_edges <- function(n, m, seed) {
  set.seed(seed)
  list(n = n, from = sample.int(n, m, replace = TRUE),
       to = sample.int(n, m, replace = TRUE), w = round(runif(m, 1, 10), 2))
}

# Out-neighbour multiset of node u as a sorted string "v:w,v:w" (order-free compare).
nbr_key_csr <- function(g, u) {
  lo <- g$offset[u]; hi <- g$offset[u + 1L] - 1L
  if (hi < lo) return(character(0))
  sort(paste(g$target[lo:hi], g$weight[lo:hi], sep = ":"))      # sort(): test side only
}
nbr_key_edges <- function(el, u) {
  idx <- el$from == u
  if (!any(idx)) return(character(0))
  sort(paste(el$to[idx], el$w[idx], sep = ":"))
}

test_that("csr_build reproduces the edge list (vs plain filtering)", {
  for (seed in 1:3) {
    el <- make_edges(30, 120, seed)
    g <- csr_build(el$n, el$from, el$to, el$w)
    expect_equal(g$m, 120L)
    expect_equal(g$offset[g$n + 1L], 121L)
    for (u in 1:el$n) expect_equal(nbr_key_csr(g, u), nbr_key_edges(el, u))
  }
})

test_that("csr_build agrees with igraph degrees and eid permutes extras", {
  el <- make_edges(50, 200, 7)
  extra <- seq_len(200)                                  # edge i carries value i
  g <- csr_build(el$n, el$from, el$to, el$w, extra = list(tag = extra))
  ig <- igraph::make_graph(rbind(el$from, el$to), n = el$n, directed = TRUE)
  deg <- igraph::degree(ig, mode = "out")
  for (u in 1:el$n) expect_equal(csr_degree(g, u), unname(deg[u]))
  # every CSR slot's tag must point back to an input edge with the same endpoints
  for (u in 1:el$n) {
    lo <- g$offset[u]; hi <- g$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) {
      i <- g$extra$tag[e]
      expect_equal(el$from[i], u); expect_equal(el$to[i], g$target[e])
    }
  }
})

test_that("csr_build edge cases: empty graph, isolated nodes, self loop", {
  g0 <- csr_build(3, integer(0), integer(0), numeric(0))
  expect_equal(g0$m, 0L); expect_equal(g0$offset, c(1L, 1L, 1L, 1L))
  g1 <- csr_build(3, 2L, 2L, 5)                          # only a self loop at node 2
  expect_equal(csr_degree(g1, 1), 0L); expect_equal(csr_degree(g1, 2), 1L)
  expect_equal(g1$target[g1$offset[2]], 2L)
})

test_that("csr_reverse swaps edge direction and keeps extras attached", {
  el <- make_edges(25, 80, 3)
  g <- csr_build(el$n, el$from, el$to, el$w, extra = list(tag = seq_len(80)))
  r <- csr_reverse(g)
  rel <- list(from = el$to, to = el$from, w = el$w)
  for (u in 1:el$n) expect_equal(nbr_key_csr(r, u), nbr_key_edges(rel, u))
  for (u in 1:el$n) {                                    # extras follow the edges
    lo <- r$offset[u]; hi <- r$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) {
      i <- r$extra$tag[e]
      expect_equal(el$to[i], u); expect_equal(el$from[i], r$target[e])
    }
  }
})

test_that("csr_bfs_subgraph returns k nodes, induced edges only", {
  el <- make_edges(60, 400, 11)
  g <- csr_build(el$n, el$from, el$to, el$w, x = runif(60), y = runif(60),
                 extra = list(tag = seq_len(400)))
  s <- csr_bfs_subgraph(g, root = 1L, k = 20L)
  orig <- attr(s, "orig_id")
  expect_lte(s$n, 20L)
  expect_equal(orig[1], 1L)
  expect_equal(s$x, g$x[orig])
  # edge count equals the number of original edges inside the chosen node set
  inside <- 0L
  for (i in seq_along(el$from)) if (el$from[i] %in% orig && el$to[i] %in% orig) inside <- inside + 1L
  expect_equal(s$m, inside)
  # asking for more nodes than exist just returns the reachable set
  big <- csr_bfs_subgraph(g, 1L, 10000L)
  expect_lte(big$n, 60L)
})

test_that("flow graph stores edges in forward/reverse pairs", {
  fg <- flow_new(4, 10)
  e1 <- flow_add_edge(fg, 1L, 2L, 7)
  e2 <- flow_add_edge(fg, 2L, 3L, 4)
  expect_equal(fg$ne, 4L)
  expect_equal(fg$to[e1], 2L); expect_equal(fg$to[rev_edge(e1)], 1L)
  expect_equal(fg$cap[e1], 7);  expect_equal(fg$cap[rev_edge(e1)], 0)
  expect_equal(rev_edge(rev_edge(e2)), e2)
  # adjacency list of node 2 holds the reverse of e1 and the forward e2
  seen <- integer(0); e <- fg$head[2]
  while (e != 0L) { seen <- c(seen, e); e <- fg$nxt[e] }
  expect_setequal(seen, c(rev_edge(e1), e2))
  # push 3 units along e1, then reset
  fg$cap[e1] <- fg$cap[e1] - 3; fg$cap[rev_edge(e1)] <- fg$cap[rev_edge(e1)] + 3
  flow_reset(fg)
  expect_equal(fg$cap[e1], 7); expect_equal(fg$cap[rev_edge(e1)], 0)
})

test_that("flow_from_csr mirrors the CSR edges", {
  el <- make_edges(15, 50, 5)
  g <- csr_build(el$n, el$from, el$to, el$w)
  fg <- flow_from_csr(g, capacity = rep(2, g$m))
  expect_equal(fg$ne, 2L * g$m)
  total_cap <- 0
  for (e in seq_len(fg$ne)) total_cap <- total_cap + fg$cap[e]
  expect_equal(total_cap, 2 * g$m)
})
