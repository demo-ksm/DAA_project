# Module 6: randomized quicksort and Karger's min cut.

# ---- quicksort -----------------------------------------------------------------
test_that("both quicksorts sort correctly on many input shapes", {
  set.seed(1)
  for (kind in c("random", "sorted", "reversed", "organ", "fewuniq")) for (n in c(1, 2, 3, 10, 257)) {
    x <- quicksort_input(n, kind)
    for (alg in list(randomized_quicksort, deterministic_quicksort)) {
      r <- alg(x)
      expect_equal(x[r$perm], sort(x), info = paste(kind, n))
      expect_equal(sort(r$perm), seq_len(n))                      # it is a permutation
    }
  }
  expect_equal(randomized_quicksort(numeric(0))$perm, integer(0))
  expect_equal(randomized_quicksort(c(3, 1, 2, 1, 3))$perm |> (\(p) c(3, 1, 2, 1, 3)[p])(), c(1, 1, 2, 3, 3))
})

test_that("deterministic quicksort is quadratic on sorted input, randomized is n log n", {
  n <- 400
  det <- deterministic_quicksort(quicksort_input(n, "sorted"))$comparisons
  expect_equal(det, n * (n - 1) / 2)                              # exactly the worst case
  set.seed(2)
  rnd <- numeric(15); for (i in 1:15) rnd[i] <- randomized_quicksort(quicksort_input(n, "sorted"))$comparisons
  expect_lt(mean(rnd), 2.2 * n * log(n))                          # ~ 2 n ln n expected
  expect_lt(mean(rnd), det / 4)
  # expected comparisons ~ 2 n ln n on random data too
  rr <- numeric(10); for (i in 1:10) rr[i] <- randomized_quicksort(runif(1000))$comparisons
  expect_gt(mean(rr), 0.8 * 2 * 1000 * (log(1000) - 1)); expect_lt(mean(rr), 1.3 * 2 * 1000 * log(1000))
})

test_that("quicksort keeps its stack logarithmic even on adversarial input", {
  det <- deterministic_quicksort(quicksort_input(2000, "sorted"))
  expect_lte(det$max_stack, ceiling(log2(2000)) + 2)
})

test_that("rank_incidents orders most urgent first", {
  sev <- c(1, 5, 3, 5, 2); dem <- c(1, 2, 9, 7, 3)
  p <- rank_incidents(sev, dem)
  expect_equal(p[1:2], c(4L, 2L))                                  # severity 5: demand 7 before 2
  expect_equal(sev[p], c(5, 5, 3, 2, 1))
  expect_equal(rank_incidents(sev, dem, randomized = FALSE), p)
})

# ---- Karger --------------------------------------------------------------------
two_cliques <- function(a, b) {                       # K_a and K_b joined by 2 bridge edges
  eu <- integer(0); ev <- integer(0)
  for (i in 1:(a - 1)) for (j in (i + 1):a) { eu <- c(eu, i); ev <- c(ev, j) }
  for (i in 1:(b - 1)) for (j in (i + 1):b) { eu <- c(eu, a + i); ev <- c(ev, a + j) }
  list(n = a + b, eu = c(eu, 1L, 2L), ev = c(ev, a + 1L, a + 2L))
}

test_that("union-find basics", {
  u <- uf_new(6)
  expect_true(uf_union(u, 1, 2)); expect_true(uf_union(u, 2, 3)); expect_false(uf_union(u, 1, 3))
  expect_equal(uf_find(u, 1), uf_find(u, 3)); expect_false(uf_find(u, 4) == uf_find(u, 1))
  expect_equal(u$count, 4L)
})

test_that("Karger finds the planted minimum cut and the cut is a real cut", {
  set.seed(6)
  g <- two_cliques(6, 6)
  r <- karger_min_cut(g$n, g$eu, g$ev, trials = 60)
  expect_equal(r$best$cut, 2L)
  expect_equal(r$best$cut, mincut_brute(g$n, g$eu, g$ev))
  s <- r$best$side
  expect_equal(sum(s != c(rep(s[1], 6), rep(!s[1], 6))), 0)           # the two cliques are separated
  for (e in r$best$cut_edges) expect_true(s[g$eu[e]] != s[g$ev[e]])
})

test_that("Karger agrees with igraph::min_cut and brute force on random graphs", {
  for (seed in 1:5) {
    set.seed(seed); n <- 10
    e <- t(utils::combn(n, 2)); e <- e[runif(nrow(e)) < 0.35, , drop = FALSE]
    ig <- igraph::make_graph(as.vector(t(e)), n = n, directed = FALSE)
    if (!igraph::is_connected(ig)) next
    ref <- igraph::min_cut(ig)
    r <- karger_min_cut(n, e[, 1], e[, 2], trials = karger_trials_needed(n, 0.001))
    expect_equal(r$best$cut, ref)
    expect_equal(mincut_brute(n, e[, 1], e[, 2]), ref)
  }
})

test_that("Karger edge cases: tree (cut 1), disconnected (cut 0), two nodes, multi-edges", {
  set.seed(7)
  expect_equal(karger_min_cut(5, c(1, 2, 3, 4), c(2, 3, 4, 5), 20)$best$cut, 1L)    # path graph
  expect_equal(karger_min_cut(4, c(1, 3), c(2, 4), 10)$best$cut, 0L)                # two components
  expect_equal(karger_min_cut(2, c(1, 1, 1), c(2, 2, 2), 3)$best$cut, 3L)           # 3 parallel edges
  expect_true(edges_connected(3, c(1, 2), c(2, 3))); expect_false(edges_connected(3, 1, 2))
})

test_that("theory helpers: success bound and number of trials", {
  expect_equal(karger_success_bound(10, 1), 2 / 90)
  expect_gt(karger_success_bound(10, karger_trials_needed(10, 0.01)), 0.99 - 1e-9)
  expect_equal(karger_trials_needed(2, 0.5), ceiling(log(2)))
})

test_that("empirical Karger success rate is at least the theoretical per-trial bound", {
  set.seed(8)
  g <- two_cliques(7, 7)                                              # n = 14, min cut 2
  n <- g$n; runs <- 400; ok <- 0
  for (i in 1:runs) if (karger_run(n, g$eu, g$ev)$cut == 2L) ok <- ok + 1
  expect_gt(ok / runs, 2 / (n * (n - 1)))
})

test_that("csr_undirected_edges collapses opposite arcs", {
  g <- csr_build(3, c(1L, 2L, 2L), c(2L, 1L, 3L), c(1, 1, 1))
  e <- csr_undirected_edges(g)
  expect_equal(length(e$eu), 2L)
})
