# =============================================================================
# Module 1 - Greedy: Fractional Knapsack
# CityShield use: distribute a DIVISIBLE emergency supply (water, food, fuel)
#                 among zones so that total benefit is maximised.
# =============================================================================

#' Fractional knapsack by greedy value-density order.
#'
#' Problem: items i have value v_i and weight w_i; a knapsack holds W. We may
#'   take any FRACTION x_i in [0,1] of an item. Maximise sum v_i x_i subject to
#'   sum w_i x_i <= W.
#' Approach (greedy): rank items by density v_i / w_i (a max-heap, no sort()),
#'   take whole items while they fit, then the largest possible fraction of the
#'   first item that does not fit.
#' Correctness (exchange argument): suppose an optimal solution takes less of a
#'   denser item a than it could while taking some of a less dense item b.
#'   Moving a small weight d from b to a changes value by d*(v_a/w_a - v_b/w_b)
#'   >= 0, so the greedy choice never hurts; repeating the swap turns any
#'   optimum into the greedy solution.
#' Time complexity: O(n log n) (n heap pushes + at most n pops).
#' Space complexity: O(n).
#'
#' @param value,weight numeric vectors (weight > 0); @param capacity total W
#' @return list(total, fraction = x_i per item, used = total weight taken)
fractional_knapsack <- function(value, weight, capacity) {
  n <- length(value)
  x <- numeric(n)
  total <- 0; used <- 0
  h <- heap_new(max(4L, n))
  for (i in seq_len(n)) if (weight[i] > 0) heap_push(h, -value[i] / weight[i], i)   # max-heap via -key
  while (!heap_empty(h) && used < capacity) {
    i <- heap_pop(h)$val
    room <- capacity - used
    if (weight[i] <= room) {                  # whole item fits
      x[i] <- 1; used <- used + weight[i]; total <- total + value[i]
    } else {                                  # take only the fraction that fits
      f <- room / weight[i]
      x[i] <- f; used <- capacity; total <- total + f * value[i]
    }
  }
  list(total = total, fraction = x, used = used)
}

#' REFERENCE: try every order of the items (greedy fill each time) and keep the
#' best. Shows that no ordering beats the density order.
#' Time complexity: O(n! * n) - tiny n only.  Space: O(n).
#' Permutations are generated iteratively (Heap's algorithm).
fractional_knapsack_brute <- function(value, weight, capacity) {
  n <- length(value)
  perm <- seq_len(n); cnt <- rep(1L, n)
  fill <- function(p) {
    used <- 0; tot <- 0
    for (i in p) {
      room <- capacity - used
      if (room <= 0) break
      f <- if (weight[i] <= room) 1 else room / weight[i]
      used <- used + f * weight[i]; tot <- tot + f * value[i]
    }
    tot
  }
  best <- fill(perm)
  i <- 1L
  while (i <= n) {                          # Heap's algorithm, iterative
    if (cnt[i] < i) {
      if (i %% 2L == 1L) { t <- perm[1L]; perm[1L] <- perm[i]; perm[i] <- t }
      else               { t <- perm[cnt[i]]; perm[cnt[i]] <- perm[i]; perm[i] <- t }
      v <- fill(perm); if (v > best) best <- v
      cnt[i] <- cnt[i] + 1L; i <- 1L
    } else { cnt[i] <- 1L; i <- i + 1L }
  }
  best
}
