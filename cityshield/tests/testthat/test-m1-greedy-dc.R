# Module 1: fractional knapsack and the divide-and-conquer maximum subarray.

# ---- fractional knapsack ------------------------------------------------------
test_that("fractional knapsack: textbook example and validity", {
  r <- fractional_knapsack(c(60, 100, 120), c(10, 20, 30), 50)
  expect_equal(r$total, 240)
  expect_equal(r$fraction, c(1, 1, 2 / 3))
  expect_equal(r$used, 50)
})

test_that("fractional knapsack matches the all-orderings brute force", {
  for (seed in 1:15) {
    set.seed(seed); n <- sample(1:6, 1)
    v <- round(runif(n, 1, 50)); w <- round(runif(n, 1, 20)); W <- round(runif(1, 1, 60))
    r <- fractional_knapsack(v, w, W)
    expect_equal(r$total, fractional_knapsack_brute(v, w, W), info = paste("seed", seed))
    expect_lte(sum(r$fraction * w), W + 1e-9)
    expect_true(all(r$fraction >= 0 & r$fraction <= 1))
    expect_equal(sum(r$fraction * v), r$total)
  }
})

test_that("fractional knapsack edge cases", {
  expect_equal(fractional_knapsack(c(5, 6), c(2, 3), 0)$total, 0)       # no capacity
  expect_equal(fractional_knapsack(numeric(0), numeric(0), 10)$total, 0)
  expect_equal(fractional_knapsack(c(5, 6), c(2, 3), 100)$total, 11)    # everything fits
  expect_equal(fractional_knapsack(c(10), c(4), 2)$total, 5)            # single item, half
  expect_equal(fractional_knapsack(c(5, 7), c(0, 2), 1)$total, 3.5)     # zero-weight item skipped
})

# ---- maximum subarray (divide and conquer) ---------------------------------------
# Test-side reference: try every (i, j) with a running sum.
msa_ref <- function(a) {
  best <- -Inf
  for (i in seq_along(a)) { s <- 0; for (j in i:length(a)) { s <- s + a[j]; if (s > best) best <- s } }
  best
}

test_that("max subarray: CLRS example and agreement with a brute-force reference", {
  a <- c(13, -3, -25, 20, -3, -16, -23, 18, 20, -7, 12, -5, -22, 15, -4, 7)
  r <- max_subarray_dc(a); expect_equal(r$sum, 43); expect_equal(c(r$start, r$end), c(8L, 11L))
  for (seed in 1:30) {
    set.seed(seed); b <- round(rnorm(sample(1:60, 1), 0, 10))
    r <- max_subarray_dc(b)
    expect_equal(r$sum, msa_ref(b))
    expect_equal(sum(b[r$start:r$end]), r$sum)                          # indices really achieve it
  }
})

test_that("max subarray edge cases: all negative, single, zeros, empty", {
  expect_equal(max_subarray_dc(c(-5, -2, -9))$sum, -2)
  expect_equal(max_subarray_dc(7)$sum, 7)
  expect_equal(max_subarray_dc(c(0, 0, 0))$sum, 0)
  expect_error(max_subarray_dc(numeric(0)))
  expect_equal(max_subarray_dc(c(1, 2, 3, 4))$sum, 10)
  expect_equal(hourly_counts(c(0, 59, 60, 125, 5000), 3), c(2, 1, 1))
})
