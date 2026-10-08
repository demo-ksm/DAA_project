# =============================================================================
# core/sort_utils.R  -  A hand-written stable sort used by geometry code
# (R's sort()/order() are banned inside graded algorithms)
# =============================================================================

#' Stable bottom-up merge sort that returns a PERMUTATION.
#'
#' Problem: order items by a numeric key without sort()/order().
#' Approach: iterative merge sort. Runs of length 1, 2, 4, ... are merged
#'   pairwise; we move an index vector, never the data. Equal keys keep their
#'   input order (stable), so sorting by a second key first and then by the
#'   main key gives a lexicographic order.
#' Time complexity: O(n log n) (log n merge rounds, each O(n)).
#' Space complexity: O(n) for the temporary index buffer.
#'
#' @param keys numeric vector (no NA)
#' @return integer vector p such that keys[p] is non-decreasing
merge_sort_perm <- function(keys) {
  n <- length(keys)
  perm <- seq_len(n)
  if (n < 2L) return(perm)
  tmp <- integer(n)
  width <- 1L
  while (width < n) {
    i <- 1L
    while (i <= n) {
      mid <- min(i + width, n + 1L)          # left run  = [i, mid)
      hi  <- min(i + 2L * width, n + 1L)     # right run = [mid, hi)
      a <- i; b <- mid; k <- i
      while (a < mid && b < hi) {
        if (keys[perm[b]] < keys[perm[a]]) { tmp[k] <- perm[b]; b <- b + 1L }
        else                               { tmp[k] <- perm[a]; a <- a + 1L }
        k <- k + 1L
      }
      while (a < mid) { tmp[k] <- perm[a]; a <- a + 1L; k <- k + 1L }
      while (b < hi)  { tmp[k] <- perm[b]; b <- b + 1L; k <- k + 1L }
      i <- i + 2L * width
    }
    t <- perm; perm <- tmp; tmp <- t         # swap buffers
    width <- width * 2L
  }
  perm
}
