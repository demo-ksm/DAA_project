# =============================================================================
# Module 5 (data structure) - Treap: the balanced tree behind the sweep line
#
# A treap is a binary search tree whose nodes also carry a random PRIORITY and
# obey the heap property on it. Random priorities make the shape that of a
# random BST, so the expected height is O(log n) without any explicit
# balancing bookkeeping (only rotations).
#
# Nodes are ENVIRONMENTS (reference semantics) with fields
#   id, pri, left, right, parent.
# The tree never compares keys itself: callers pass a comparison closure, so the
# same tree can order segments by "y at the current sweep position".
#
# Operation costs (expected): insert, delete, succ, pred, lower_bound: O(log n)
# Space: O(n).
# =============================================================================

#' @param max_id largest id that will ever be inserted (ids are 1..max_id)
treap_new <- function(max_id) {
  tr <- new.env(parent = emptyenv())
  tr$root <- NULL
  tr$node_of <- vector("list", max_id)      # id -> node environment
  tr$size <- 0L
  tr
}

.treap_rotate_up <- function(tr, x) {
  p <- x$parent
  g <- p$parent
  if (identical(p$left, x)) {               # right rotation
    b <- x$right
    p$left <- b
    if (!is.null(b)) b$parent <- p
    x$right <- p
  } else {                                  # left rotation
    b <- x$left
    p$right <- b
    if (!is.null(b)) b$parent <- p
    x$left <- p
  }
  p$parent <- x
  x$parent <- g
  if (is.null(g)) tr$root <- x
  else if (identical(g$left, p)) g$left <- x
  else g$right <- x
}

#' Insert id. `goes_left(other_id)` must return TRUE when the NEW element sorts
#' before `other_id`. Expected O(log n).
treap_insert <- function(tr, id, goes_left) {
  node <- new.env(parent = emptyenv())
  node$id <- id; node$pri <- runif(1)
  node$left <- NULL; node$right <- NULL; node$parent <- NULL
  tr$node_of[[id]] <- node
  tr$size <- tr$size + 1L
  if (is.null(tr$root)) { tr$root <- node; return(invisible(tr)) }
  cur <- tr$root
  repeat {
    if (goes_left(cur$id)) {
      if (is.null(cur$left)) { cur$left <- node; break }
      cur <- cur$left
    } else {
      if (is.null(cur$right)) { cur$right <- node; break }
      cur <- cur$right
    }
  }
  node$parent <- cur
  while (!is.null(node$parent) && node$pri < node$parent$pri) .treap_rotate_up(tr, node)
  invisible(tr)
}

#' Delete element `id` (located through its node, no key comparison needed):
#' rotate it down until it is a leaf, then cut it off. Expected O(log n).
treap_delete <- function(tr, id) {
  x <- tr$node_of[[id]]
  if (is.null(x)) stop("id not in treap")
  repeat {
    l <- x$left; r <- x$right
    if (is.null(l) && is.null(r)) break
    c <- if (is.null(l)) r else if (is.null(r)) l else if (l$pri < r$pri) l else r
    .treap_rotate_up(tr, c)                 # lifts c above x, pushing x down
  }
  p <- x$parent
  if (is.null(p)) tr$root <- NULL
  else if (identical(p$left, x)) p$left <- NULL
  else p$right <- NULL
  tr$node_of[id] <- list(NULL)
  tr$size <- tr$size - 1L
  invisible(tr)
}

#' In-order successor id (0 if none). O(log n) expected.
treap_succ <- function(tr, id) {
  x <- tr$node_of[[id]]
  if (!is.null(x$right)) {
    y <- x$right
    while (!is.null(y$left)) y <- y$left
    return(y$id)
  }
  y <- x
  while (!is.null(y$parent) && identical(y$parent$right, y)) y <- y$parent
  if (is.null(y$parent)) 0L else y$parent$id
}

#' In-order predecessor id (0 if none). O(log n) expected.
treap_pred <- function(tr, id) {
  x <- tr$node_of[[id]]
  if (!is.null(x$left)) {
    y <- x$left
    while (!is.null(y$right)) y <- y$right
    return(y$id)
  }
  y <- x
  while (!is.null(y$parent) && identical(y$parent$left, y)) y <- y$parent
  if (is.null(y$parent)) 0L else y$parent$id
}

#' Largest element id (0 if empty). O(log n) expected.
treap_max <- function(tr) {
  y <- tr$root
  if (is.null(y)) return(0L)
  while (!is.null(y$right)) y <- y$right
  y$id
}

#' First element (in order) for which `ge(id)` is TRUE; 0 if none.
#' `ge` must be monotone along the in-order sequence (FALSE...FALSE TRUE...TRUE).
#' O(log n) expected.
treap_lower_bound <- function(tr, ge) {
  y <- tr$root
  best <- 0L
  while (!is.null(y)) {
    if (ge(y$id)) { best <- y$id; y <- y$left } else y <- y$right
  }
  best
}

#' In-order list of ids using an explicit stack (no recursion). O(n).
treap_inorder <- function(tr) {
  out <- integer(tr$size); k <- 0L
  st <- list(); top <- 0L
  y <- tr$root
  while (!is.null(y) || top > 0L) {
    while (!is.null(y)) { top <- top + 1L; st[[top]] <- y; y <- y$left }
    y <- st[[top]]; top <- top - 1L
    k <- k + 1L; out[k] <- y$id
    y <- y$right
  }
  out
}

#' Height of the tree (levels), iterative BFS. O(n). Used in tests/benchmarks to
#' show the treap stays logarithmic.
treap_height <- function(tr) {
  if (is.null(tr$root)) return(0L)
  level <- list(tr$root); h <- 0L
  while (length(level) > 0L) {
    h <- h + 1L
    nxt <- list()
    for (nd in level) {
      if (!is.null(nd$left))  nxt[[length(nxt) + 1L]] <- nd$left
      if (!is.null(nd$right)) nxt[[length(nxt) + 1L]] <- nd$right
    }
    level <- nxt
  }
  h
}
