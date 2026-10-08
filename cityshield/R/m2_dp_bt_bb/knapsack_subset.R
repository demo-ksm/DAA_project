# =============================================================================
# Module 2 - DP: 0/1 Knapsack
# CityShield use: pick infrastructure UPGRADES (each all-or-nothing) under a
#                 budget to maximise resilience benefit.
# =============================================================================

#' 0/1 knapsack by dynamic programming.
#'
#' Problem: items with integer weight w_i and value v_i; capacity W; take each
#'   item at most once; maximise total value.
#' Approach (DP): K[i,c] = best value using the first i items with capacity c.
#'   K[i,c] = max(K[i-1,c], K[i-1,c-w_i] + v_i). Either item i is skipped or it
#'   is taken, and the rest is an optimal solution of a smaller instance.
#'   Walk back through the table to find which items were taken.
#' Time complexity: Theta(n W) - PSEUDO-polynomial (exponential in the number of
#'   bits of W).
#' Space complexity: Theta(n W) for the table.
#' @return list(value, items = indices taken, weight)
knapsack01 <- function(w, v, W) {
  n <- length(w)
  K <- matrix(0, n + 1L, W + 1L)
  for (i in seq_len(n)) {
    for (c in 0:W) {
      best <- K[i, c + 1L]
      if (w[i] <= c) { alt <- K[i, c - w[i] + 1L] + v[i]; if (alt > best) best <- alt }
      K[i + 1L, c + 1L] <- best
    }
  }
  items <- integer(0); c <- W
  for (i in rev(seq_len(n))) {
    if (K[i + 1L, c + 1L] != K[i, c + 1L]) { items <- c(i, items); c <- c - w[i] }
  }
  list(value = K[n + 1L, W + 1L], items = items, weight = sum_idx(w, items))
}

#' REFERENCE: try all 2^n subsets. Time O(n 2^n), space O(n).
knapsack01_brute <- function(w, v, W) {
  n <- length(w); best <- 0
  for (mask in 0:(2^n - 1)) {
    x <- mask; tw <- 0; tv <- 0
    for (i in seq_len(n)) { if (x %% 2 == 1) { tw <- tw + w[i]; tv <- tv + v[i] }; x <- x %/% 2 }
    if (tw <= W && tv > best) best <- tv
  }
  best
}

sum_idx <- function(x, idx) { s <- 0; for (i in idx) s <- s + x[i]; s }
