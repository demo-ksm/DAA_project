# Module 2 (branch and bound): the generic engine and the TSP solver built on it.
# (csr_undirected_adj, a graph helper used by Karger and the sensor placement, is checked here too.)

test_that("csr_undirected_adj drops self loops and duplicate links", {
  g <- csr_build(3, c(1L, 2L, 1L, 3L), c(2L, 1L, 3L, 3L), c(1, 1, 1, 1))   # 1<->2 duplicate, 1->3, self loop 3->3
  a <- csr_undirected_adj(g)
  expect_setequal(a[[1]], c(2L, 3L)); expect_equal(a[[2]], 1L); expect_equal(a[[3]], 1L)
})

# ---- B&B engine + TSP ---------------------------------------------------------
rand_metric <- function(n, seed) {
  set.seed(seed); P <- cbind(runif(n, 0, 100), runif(n, 0, 100)); D <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) D[i, j] <- sqrt(sum((P[i, ] - P[j, ])^2))
  D
}

test_that("TSP B&B (all three strategies) = brute force", {
  for (seed in 1:6) {
    n <- 4 + seed %% 4; D <- rand_metric(n, seed)
    ref <- tsp_brute(D)$cost
    for (s in c("lifo", "fifo", "lc")) {
      r <- tsp_bb(D, strategy = s)
      expect_equal(r$cost, ref, info = paste(s, seed))
      expect_equal(sort(r$tour), 1:n)
      expect_equal(tour_cost(D, r$tour), r$cost)
    }
  }
})

test_that("TSP B&B on asymmetric costs, tiny cases, and node counts", {
  set.seed(3); n <- 7; D <- matrix(sample(1:50, n * n, TRUE), n, n); diag(D) <- 0
  expect_equal(tsp_bb(D)$cost, tsp_brute(D)$cost)
  expect_equal(tsp_bb(matrix(0, 1, 1))$cost, 0)
  D2 <- matrix(c(0, 3, 3, 0), 2, 2); expect_equal(tsp_bb(D2)$cost, 6)
  D8 <- rand_metric(8, 11)
  lc <- tsp_bb(D8, "lc"); lifo <- tsp_bb(D8, "lifo")
  expect_lt(lc$expanded, factorial(7))                                   # pruning beats full enumeration
  expect_lte(lc$expanded, lifo$expanded)                                 # best-first never expands more here
  cap <- tsp_bb(D8, "lifo", max_nodes = 3)
  expect_false(cap$complete)
})

test_that("tsp_reduce lowers every row/column minimum to zero", {
  M <- matrix(c(Inf, 10, 20, 30, Inf, 5, 15, 25, Inf), 3, 3, byrow = TRUE)
  r <- tsp_reduce(M)
  for (i in 1:3) expect_equal(min(r$M[i, ]), 0)
  for (j in 1:3) expect_equal(min(r$M[, j]), 0)
  expect_equal(r$red, 10 + 5 + 15)
})

test_that("generic branch_and_bound engine: strategies agree on a toy problem", {
  # choose 3 numbers (one per level) minimising their sum; bound = sum so far + min remaining
  vals <- list(c(5, 2, 9), c(4, 8, 1), c(7, 3, 6))
  root <- list(lvl = 0L, s = 0)
  expand <- function(nd) { k <- list(); for (x in vals[[nd$lvl + 1L]]) k[[length(k) + 1L]] <- list(lvl = nd$lvl + 1L, s = nd$s + x); k }
  bound <- function(nd) { b <- nd$s; if (nd$lvl < 3) for (l in (nd$lvl + 1L):3) b <- b + min(vals[[l]]); b }
  for (s in c("lifo", "fifo", "lc")) {
    r <- branch_and_bound(root, expand, bound, function(nd) nd$lvl == 3, function(nd) nd$s, strategy = s)
    expect_equal(r$cost, 2 + 1 + 3)
  }
  expect_error(branch_and_bound(root, expand, bound, function(nd) TRUE, function(nd) 0, strategy = "xyz"))
})

test_that("TSP B&B with an initial upper bound from the 2-approximation still finds the optimum", {
  for (seed in 1:5) {
    D <- rand_metric(7, 20 + seed)
    ap <- tsp_approx(D)
    r <- tsp_bb(D, "lc", ub = ap$cost + 1e-6)
    expect_equal(r$cost, tsp_brute(D)$cost)
    expect_true(r$complete)
  }
})
