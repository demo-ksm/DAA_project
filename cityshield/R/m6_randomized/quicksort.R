# =============================================================================
# Module 6 - Randomized algorithms: Quicksort
# CityShield use: rank incidents by priority (severity, then demand) so the
#   control room always sees the most urgent first - even if the incident feed
#   arrives already sorted, which is the worst case for naive quicksort.
# Both versions are ITERATIVE (explicit stack, smaller side first), so the
# stack never exceeds O(log n) frames even on adversarial inputs.
# =============================================================================

#' Quicksort a numeric key vector (returns the sorting permutation).
#'
#' Problem: order n keys ascending; we return the permutation `perm` with
#'   keys[perm] non-decreasing, plus the number of key comparisons made.
#' Approach: Lomuto partition around a PIVOT, then recurse on both sides.
#'   * deterministic: pivot = LAST element of the range. On an already sorted
#'     (or reverse sorted) input every partition is maximally unbalanced.
#'   * randomized: pivot = a uniformly random element of the range (swapped to
#'     the end first), so no fixed input is bad: the pivot's rank is uniform.
#' Expected comparisons (randomized, ANY input): 2 n ln n ~ 1.39 n log2 n.
#'   Proof sketch: elements i < j (in sorted order) are compared iff one of them
#'   is the first pivot chosen among z_i..z_j, probability 2/(j-i+1);
#'   summing over all pairs gives 2 n H_n = O(n log n).
#' Worst case: Theta(n^2) comparisons (deterministic on sorted input; for the
#'   randomized version only with probability ~ 1/n!).
#' Space complexity: O(log n) stack (always the smaller part is processed first
#'   from an explicit stack) + O(n) for the permutation.
#'
#' @param randomized TRUE = random pivot, FALSE = last-element pivot
#' @return list(perm, comparisons, max_stack)
quicksort_perm <- function(keys, randomized = TRUE) {
  n <- length(keys)
  perm <- seq_len(n)
  comps <- 0
  if (n < 2L) return(list(perm = perm, comparisons = 0, max_stack = 0L))
  stk_lo <- integer(64L); stk_hi <- integer(64L)
  top <- 1L; stk_lo[1L] <- 1L; stk_hi[1L] <- n
  max_stack <- 1L
  while (top > 0L) {
    lo <- stk_lo[top]; hi <- stk_hi[top]; top <- top - 1L
    while (lo < hi) {
      if (randomized) {                       # move a random pivot to the end
        r <- lo + sample.int(hi - lo + 1L, 1L) - 1L
        t <- perm[r]; perm[r] <- perm[hi]; perm[hi] <- t
      }
      pv <- keys[perm[hi]]
      i <- lo - 1L
      for (j in lo:(hi - 1L)) {               # Lomuto partition
        comps <- comps + 1
        if (keys[perm[j]] <= pv) { i <- i + 1L; t <- perm[i]; perm[i] <- perm[j]; perm[j] <- t }
      }
      t <- perm[i + 1L]; perm[i + 1L] <- perm[hi]; perm[hi] <- t
      p <- i + 1L
      # push the LARGER side, continue with the smaller one (bounded stack)
      if (p - lo < hi - p) {
        if (p + 1L < hi) {
          top <- top + 1L
          if (top > length(stk_lo)) { length(stk_lo) <- 2L * top; length(stk_hi) <- 2L * top }
          stk_lo[top] <- p + 1L; stk_hi[top] <- hi
        }
        hi <- p - 1L
      } else {
        if (lo < p - 1L) {
          top <- top + 1L
          if (top > length(stk_lo)) { length(stk_lo) <- 2L * top; length(stk_hi) <- 2L * top }
          stk_lo[top] <- lo; stk_hi[top] <- p - 1L
        }
        lo <- p + 1L
      }
      if (top > max_stack) max_stack <- top
    }
  }
  list(perm = perm, comparisons = comps, max_stack = max_stack)
}

#' Randomized quicksort (random pivot). See quicksort_perm().
randomized_quicksort <- function(keys) quicksort_perm(keys, randomized = TRUE)

#' Deterministic quicksort (last-element pivot). See quicksort_perm().
deterministic_quicksort <- function(keys) quicksort_perm(keys, randomized = FALSE)

#' Rank incidents most-urgent-first: key = severity * 1000 + demand (higher =
#' more urgent), sorted DEscending by quicksorting the negated key. O(n log n)
#' expected. Returns incident row indices in priority order.
rank_incidents <- function(severity, demand, randomized = TRUE) {
  key <- numeric(length(severity))
  for (i in seq_along(severity)) key[i] <- -(severity[i] * 1000 + demand[i])
  quicksort_perm(key, randomized)$perm
}

#' Adversarial / structured inputs for the benchmark. n values.
quicksort_input <- function(n, kind = "random") {
  switch(kind,
         random   = runif(n),
         sorted   = as.numeric(seq_len(n)),
         reversed = as.numeric(n:1),
         organ    = c(as.numeric(seq_len(n %/% 2)), as.numeric((n - n %/% 2):1)),   # up then down
         fewuniq  = sample.int(5L, n, replace = TRUE) + runif(n) * 1e-6,
         stop("unknown kind"))
}
