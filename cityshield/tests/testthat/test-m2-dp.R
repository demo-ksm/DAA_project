# Module 2 (dynamic programming part): longest common subsequence and 0/1 knapsack.

# ---- LCS ---------------------------------------------------------------------
test_that("LCS: classic examples, validity, and brute force", {
  a <- strsplit("ABCBDAB", "")[[1]]; b <- strsplit("BDCABA", "")[[1]]
  r <- lcs(a, b); expect_equal(r$length, 4L)
  is_subseq <- function(s, x) { p <- 1L; for (ch in x) if (p <= length(s) && ch == s[p]) p <- p + 1L; p > length(s) }
  expect_true(is_subseq(r$subsequence, a)); expect_true(is_subseq(r$subsequence, b))
  expect_equal(length(r$subsequence), 4L)
  expect_equal(lcs_length(a, b), 4L)
  set.seed(2)
  for (i in 1:25) {
    x <- sample(letters[1:4], sample(0:10, 1), TRUE); y <- sample(letters[1:4], sample(0:10, 1), TRUE)
    ref <- lcs_brute(x, y)
    expect_equal(lcs(x, y)$length, ref); expect_equal(lcs_length(x, y), ref)
    expect_equal(lcs_length(y, x), ref)                               # symmetric
  }
})

test_that("LCS edge cases: empty, identical, disjoint", {
  expect_equal(lcs(character(0), c("a"))$length, 0L)
  expect_equal(lcs_length(character(0), character(0)), 0L)
  expect_equal(lcs(c("a", "b"), c("a", "b"))$length, 2L)
  expect_equal(lcs(c("a", "b"), c("c", "d"))$length, 0L)
  expect_equal(lcs(1:5, 3:8)$subsequence, 3:5)                        # integers work too
})

test_that("report similarity and de-duplication find the planted duplicates", {
  expect_equal(report_similarity("fire at Adyar road", "FIRE at adyar ROAD"), 1)
  expect_lt(report_similarity("fire at adyar", "flood in tambaram"), 0.3)
  skip_if_not(file.exists(data_path("synthetic", "call_logs.rds")))
  cl <- readRDS(data_path("synthetic", "call_logs.rds"))
  sub <- cl[c(1:120, which(cl$is_duplicate & cl$dup_of <= 120)), ]
  b <- report_best_matches(sub$text, 0.5)
  truth <- sub$is_duplicate
  # a trade-off, not a free lunch: low thresholds flag template-similar reports, high ones miss noisy copies
  pr <- b$sim >= 0.9
  recall <- sum(truth & pr) / sum(truth); precision <- sum(truth & pr) / max(1, sum(pr))
  expect_gt(recall, 0.85); expect_gt(precision, 0.7)
  lo <- b$sim >= 0.7
  expect_gt(sum(truth & lo) / sum(truth), recall - 1e-9)             # lower threshold never loses recall
  expect_lt(sum(truth & lo) / max(1, sum(lo)), precision)            # ... but hurts precision
  # the planted duplicate really points at its original (best match = the source report)
  hit <- 0; tot <- 0
  for (k in which(truth)) { tot <- tot + 1; if (b$best_i[k] == sub$dup_of[k]) hit <- hit + 1 }
  expect_gt(hit / tot, 0.9)
  d <- dedupe_reports(sub$text[1:60], threshold = 0.9)
  expect_true(all(is.na(d) | d < seq_along(d)))                       # duplicates always point to EARLIER reports
})

# ---- 0/1 knapsack ------------------------------------------------------------
test_that("0/1 knapsack: textbook value, reconstruction and brute force", {
  r <- knapsack01(c(10, 20, 30), c(60, 100, 120), 50)
  expect_equal(r$value, 220); expect_equal(r$items, c(2L, 3L))
  for (seed in 1:25) {
    set.seed(seed); n <- sample(1:10, 1)
    w <- sample(1:15, n, TRUE); v <- sample(1:40, n, TRUE); W <- sample(5:40, 1)
    r <- knapsack01(w, v, W)
    expect_equal(r$value, knapsack01_brute(w, v, W), info = paste(seed))
    expect_lte(sum(w[r$items]), W); expect_equal(sum(v[r$items]), r$value)
  }
})

test_that("0/1 knapsack edge cases", {
  expect_equal(knapsack01(c(5), c(10), 4)$value, 0)                    # nothing fits
  expect_equal(knapsack01(c(5), c(10), 5)$items, 1L)
  expect_equal(knapsack01(integer(0), numeric(0), 10)$value, 0)
  expect_equal(knapsack01(c(1, 1), c(3, 4), 0)$value, 0)
})
