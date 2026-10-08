# =============================================================================
# Module 4 - Maximum bipartite matching via max flow
# CityShield use: assign ambulances (left) to incidents (right) so that as many
#                 incidents as possible get exactly one reachable ambulance.
# =============================================================================

#' Maximum bipartite matching through a unit-capacity flow network.
#'
#' Problem: given left nodes 1..nl, right nodes 1..nr and "compatible" pairs
#'   (l_i, r_i), select the largest set of pairs that share no endpoint.
#' Approach: reduce to max flow. Network: source -> every left node (cap 1),
#'   left -> right for each compatible pair (cap 1), every right node -> sink
#'   (cap 1). An integral max flow of value k is a matching of size k
#'   (flow integrality theorem); matched pairs = middle edges carrying 1.
#'   The flow is computed with push-relabel (run to completion so that the
#'   edge flows are a valid flow and the pairs can be read off).
#' Time complexity: O(V^3) for FIFO push-relabel (V = nl + nr + 2); much faster
#'   in practice on these sparse unit-capacity networks.
#' Space complexity: O(nl + nr + E).
#'
#' @param nl,nr number of left / right nodes
#' @param l,r   integer vectors, pair i is (l[i], r[i])
#' @return list(size, match_left = right partner of each left node (0 = none),
#'              match_right, pairs = 2-column matrix)
max_bipartite_matching <- function(nl, nr, l, r) {
  m <- length(l)
  s <- nl + nr + 1L; t <- nl + nr + 2L
  fg <- flow_new(nl + nr + 2L, m + nl + nr)
  for (i in seq_len(nl)) flow_add_edge(fg, s, i, 1)
  mid <- integer(m)
  for (k in seq_len(m)) mid[k] <- flow_add_edge(fg, l[k], nl + r[k], 1)
  for (j in seq_len(nr)) flow_add_edge(fg, nl + j, t, 1)

  res <- push_relabel(fg, s, t, complete = TRUE)    # need a VALID flow to read the pairs

  match_left <- integer(nl); match_right <- integer(nr)
  for (k in seq_len(m)) {
    if (fg$cap[mid[k]] == 0) {               # saturated middle edge = matched pair
      match_left[l[k]] <- r[k]; match_right[r[k]] <- l[k]
    }
  }
  size <- 0L
  pairs <- matrix(0L, nrow = min(nl, nr), ncol = 2L)
  for (i in seq_len(nl)) if (match_left[i] != 0L) { size <- size + 1L; pairs[size, ] <- c(i, match_left[i]) }
  list(size = size, match_left = match_left, match_right = match_right,
       pairs = pairs[seq_len(size), , drop = FALSE], flow = res$value)
}

#' Brute-force maximum matching (REFERENCE) by exhaustive backtracking over the
#' left nodes with an EXPLICIT stack (no recursion). For tiny instances only.
#'
#' Problem: same as max_bipartite_matching().
#' Approach: for each left node try every unused compatible right node, or
#'   skip it; remember the best total.
#' Time complexity: O((d+1)^nl) where d = max degree.  Space: O(nl).
matching_brute <- function(nl, nr, l, r) {
  adj <- vector("list", nl)
  for (k in seq_along(l)) adj[[l[k]]] <- c(adj[[l[k]]], r[k])
  used <- logical(nr)
  choice <- integer(nl)                      # index into adj[[i]]; deg+1 = skip
  best <- 0L; size <- 0L
  i <- 1L; choice[1L] <- 0L
  if (nl == 0L) return(0L)
  repeat {
    # advance choice at level i
    choice[i] <- choice[i] + 1L
    d <- length(adj[[i]])
    if (choice[i] > 1L + d) {                 # level exhausted: backtrack
      i <- i - 1L
      if (i == 0L) break
      # undo the choice at level i (if it matched something)
      if (choice[i] >= 1L && choice[i] <= length(adj[[i]])) {
        used[adj[[i]][choice[i]]] <- FALSE; size <- size - 1L
      }
      next
    }
    if (choice[i] <= d) {
      rr <- adj[[i]][choice[i]]
      if (used[rr]) next                      # right node already taken
      used[rr] <- TRUE; size <- size + 1L
    }                                         # choice == d+1 means "skip i"
    if (size > best) best <- size
    if (i == nl) {                            # leaf: undo own choice and retry
      if (choice[i] <= d) { used[adj[[i]][choice[i]]] <- FALSE; size <- size - 1L }
      next
    }
    i <- i + 1L; choice[i] <- 0L
  }
  best
}
