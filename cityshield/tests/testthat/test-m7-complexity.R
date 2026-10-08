# Module 7: vertex cover, set cover and the TSP approximation (exact vs approximate).

# ---- vertex cover -------------------------------------------------------------------
test_that("vertex cover: approx is valid and <= 2 OPT; exact = brute force = n - alpha (igraph)", {
  for (seed in 1:15) {
    set.seed(seed); n <- sample(5:14, 1)
    e <- t(utils::combn(n, 2)); e <- e[runif(nrow(e)) < 0.3, , drop = FALSE]
    if (nrow(e) == 0) next
    ap <- vertex_cover_approx(n, e[, 1], e[, 2])
    ex <- vertex_cover_exact(n, e[, 1], e[, 2])
    br <- vertex_cover_brute(n, e[, 1], e[, 2])
    expect_true(is_vertex_cover(e[, 1], e[, 2], ap$cover, n)); expect_true(is_vertex_cover(e[, 1], e[, 2], ex$cover, n))
    expect_equal(ex$size, br, info = paste(seed)); expect_true(ex$complete)
    expect_lte(ap$size, 2 * br); expect_gte(ap$size, br)
    ig <- igraph::make_graph(as.vector(t(e)), n = n, directed = FALSE)
    expect_equal(br, n - igraph::independence_number(ig))             # VC and IS are complements
    expect_equal(ap$size, 2 * ap$matching)
  }
})

test_that("vertex cover edge cases", {
  expect_equal(vertex_cover_exact(3, integer(0), integer(0))$size, 0L)
  expect_equal(vertex_cover_approx(3, integer(0), integer(0))$size, 0L)
  star <- vertex_cover_exact(6, rep(1L, 5), 2:6); expect_equal(star$size, 1L)
  expect_equal(vertex_cover_approx(6, rep(1L, 5), 2:6)$size, 2L)     # approx is not optimal on a star (ratio exactly 2)
  tri <- vertex_cover_exact(3, c(1L, 2L, 1L), c(2L, 3L, 3L)); expect_equal(tri$size, 2L)
  expect_false(vertex_cover_exact(12, c(1:11), c(2:12), max_nodes = 0)$complete)
})

# ---- set cover ------------------------------------------------------------------------
test_that("set cover: exact = brute force; greedy valid and within H(d) * OPT", {
  for (seed in 1:15) {
    set.seed(seed); U <- sample(6:14, 1); S <- sample(5:11, 1)
    sets <- vector("list", S)
    for (s in 1:S) sets[[s]] <- sort(sample.int(U, sample(2:5, 1)))
    cov <- logical(U); for (s in sets) cov[s] <- TRUE
    if (!all(cov)) sets[[S + 1L]] <- which(!cov)                      # make it feasible
    gr <- set_cover_greedy(sets, U); ex <- set_cover_exact(sets, U); br <- set_cover_brute(sets, U)
    expect_true(gr$covered_all); expect_equal(ex$size, br, info = paste(seed))
    d <- max(lengths(sets))
    expect_gte(gr$size, br); expect_lte(gr$size, harmonic_bound(d) * br + 1e-9)
    covered <- logical(U); for (s in gr$chosen) covered[sets[[s]]] <- TRUE
    expect_true(all(covered))
    covered <- logical(U); for (s in ex$chosen) covered[sets[[s]]] <- TRUE
    expect_true(all(covered))
  }
})

test_that("set cover: greedy can be suboptimal (classic instance) and infeasible detection", {
  # universe 1..6: greedy picks the 3-set first and then needs 2 more; OPT = 2
  sets <- list(c(1L, 2L, 3L, 4L), c(1L, 2L, 5L), c(3L, 4L, 6L), c(5L, 6L))
  expect_equal(set_cover_exact(sets, 6L)$size, 2L)
  expect_gte(set_cover_greedy(sets, 6L)$size, 2L)
  expect_false(set_cover_greedy(list(1L), 2L)$covered_all)
  expect_true(is.na(set_cover_exact(list(1L), 2L)$size))
  expect_equal(set_cover_greedy(list(), 0L)$size, 0L)
  expect_equal(harmonic_bound(4), 1 + 1/2 + 1/3 + 1/4)
  D <- matrix(c(0, 1, 5, 1, 0, 1, 5, 1, 0), 3, 3)
  expect_equal(coverage_sets(D, 1), list(1:2, 1:3, 2:3))
})

# ---- TSP approximation ----------------------------------------------------------------
test_that("MST (Prim) weight matches igraph and the 2-approx tour is valid and within 2x optimal", {
  for (seed in 1:10) {
    set.seed(seed); n <- sample(5:8, 1)
    P <- cbind(runif(n), runif(n)); D <- as.matrix(dist(P))
    mst <- mst_prim(D)
    ref <- sum(igraph::E(igraph::mst(igraph::graph_from_adjacency_matrix(D, mode = "undirected", weighted = TRUE)))$weight)
    expect_equal(mst$weight, ref)
    ap <- tsp_approx(D)
    expect_equal(sort(ap$tour), 1:n)
    expect_equal(ap$cost, tour_cost(D, ap$tour))
    opt <- tsp_brute(D)$cost
    expect_gte(ap$cost, opt - 1e-9); expect_lte(ap$cost, 2 * opt + 1e-9)
    expect_lte(mst$weight, opt + 1e-9)                                  # MST <= OPT
    expect_lte(ap$cost, 2 * mst$weight + 1e-9)
  }
  expect_equal(tsp_approx(matrix(0, 1, 1))$cost, 0)
  expect_equal(tsp_approx(matrix(c(0, 4, 4, 0), 2, 2))$cost, 8)
})
