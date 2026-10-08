# =============================================================================
# Module 1 - Divide and Conquer: Maximum subarray
# CityShield use: find the worst CONTINUOUS period of incident load - the
#                 stretch of consecutive hours whose (load - baseline capacity)
#                 sum is largest.
# =============================================================================

#' Maximum subarray by divide and conquer.
#'
#' Problem: in a numeric vector A find a contiguous block with the largest sum.
#' Approach: split at the middle. The best block is (a) entirely in the left
#'   half, (b) entirely in the right half, or (c) CROSSING the middle. Case (c)
#'   is found directly: best suffix of the left half + best prefix of the right
#'   half (two linear scans). Recurse on (a) and (b).
#'   Recursion depth is only log2(n), so plain recursion is safe here.
#' Recurrence: T(n) = 2T(n/2) + Theta(n)  =>  Theta(n log n) (Master theorem).
#' Space complexity: O(log n) recursion stack.
#' Correctness: the three cases are exhaustive and (c) is optimal because the
#'   crossing block splits into an independent best suffix + best prefix.
#'
#' @return list(sum, start, end)
max_subarray_dc <- function(a) {
  rec <- function(lo, hi) {
    if (lo == hi) return(list(sum = a[lo], start = lo, end = lo))
    mid <- (lo + hi) %/% 2L
    L <- rec(lo, mid); R <- rec(mid + 1L, hi)
    # best block crossing the midpoint
    s <- 0; best_l <- -Inf; bl <- mid
    for (i in mid:lo) { s <- s + a[i]; if (s > best_l) { best_l <- s; bl <- i } }
    s <- 0; best_r <- -Inf; br <- mid + 1L
    for (j in (mid + 1L):hi) { s <- s + a[j]; if (s > best_r) { best_r <- s; br <- j } }
    C <- list(sum = best_l + best_r, start = bl, end = br)
    if (L$sum >= R$sum && L$sum >= C$sum) L else if (R$sum >= C$sum) R else C
  }
  if (length(a) == 0L) stop("empty input")
  rec(1L, length(a))
}

#' Hourly incident load per hour (counts from incident time stamps in minutes).
#' O(n + H). Used to build the max-subarray input: load - capacity.
hourly_counts <- function(time_min, n_hours) {
  out <- numeric(n_hours)
  for (t in time_min) {
    h <- t %/% 60 + 1
    if (h >= 1 && h <= n_hours) out[h] <- out[h] + 1
  }
  out
}
