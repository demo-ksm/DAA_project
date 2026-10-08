# Validation of Dijkstra (and the all-pairs metric built from it) against igraph.

test_that("dijkstra matches igraph and the O(V^2) array version", {
  for (seed in 1:4) {
    el <- rand_edges(40, 160, seed)
    g <- csr_build(el$n, el$from, el$to, el$w)
    ig <- to_igraph(el)
    for (s in c(1L, 7L, 40L)) {
      d <- dijkstra(g, s)$dist
      ref <- as.numeric(igraph::distances(ig, v = s, mode = "out", algorithm = "dijkstra"))
      expect_equal(d, ref)
      expect_equal(dijkstra_array(g, s)$dist, ref)
    }
  }
})

test_that("dijkstra: multi-source, target early exit, limit, unreachable", {
  el <- rand_edges(30, 90, 5)
  g <- csr_build(el$n, el$from, el$to, el$w)
  d1 <- dijkstra(g, 3L)$dist; d2 <- dijkstra(g, 9L)$dist
  dm <- dijkstra(g, c(3L, 9L))$dist
  expect_equal(dm, pmin(d1, d2))                       # pmin(): test-side check only
  full <- dijkstra(g, 3L)
  tgt <- 20L
  early <- dijkstra(g, 3L, target = tgt)
  if (is.finite(full$dist[tgt])) expect_equal(early$dist[tgt], full$dist[tgt])
  expect_lte(early$settled, full$settled)
  lim <- dijkstra(g, 3L, limit = 8)
  inside <- full$dist <= 8
  expect_equal(lim$dist[inside], full$dist[inside])
  # isolated node
  g2 <- csr_build(3, 1L, 2L, 1)
  expect_equal(dijkstra(g2, 1L)$dist, c(0, 1, Inf))
  expect_equal(dijkstra(g2, 3L)$dist, c(Inf, Inf, 0))
})

test_that("path_from_parent returns a path whose weight equals the distance", {
  el <- rand_edges(30, 120, 8)
  g <- csr_build(el$n, el$from, el$to, el$w)
  r <- dijkstra(g, 1L)
  for (dst in 2:30) {
    p <- path_from_parent(r$parent, 1L, dst, r$dist)
    if (!is.finite(r$dist[dst])) { expect_length(p, 0); next }
    expect_equal(p[1], 1L); expect_equal(p[length(p)], dst)
    total <- 0
    for (i in seq_len(length(p) - 1L)) {                # cheapest parallel edge p[i] -> p[i+1]
      best <- Inf
      lo <- g$offset[p[i]]
      for (k in seq_len(g$offset[p[i] + 1L] - lo))
        if (g$target[lo + k - 1L] == p[i + 1L] && g$weight[lo + k - 1L] < best) best <- g$weight[lo + k - 1L]
      total <- total + best
    }
    expect_equal(total, r$dist[dst])
  }
})

test_that("dijkstra single node / no edges", {
  g <- csr_build(1, integer(0), integer(0), numeric(0))
  expect_equal(dijkstra(g, 1L)$dist, 0)
})

test_that("metric_closure = symmetrised all-pairs Dijkstra distances (checked against igraph)", {
  el <- rand_edges(25, 120, 12)
  g <- csr_build(el$n, el$from, el$to, el$w)
  S <- metric_closure(g)
  ref <- igraph::distances(to_igraph(el), mode = "out", algorithm = "dijkstra")
  for (i in 1:25) for (j in 1:25) {
    expect_equal(S[i, j], (ref[i, j] + ref[j, i]) / 2)
    expect_equal(S[i, j], S[j, i])
  }
})
