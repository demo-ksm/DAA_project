# =============================================================================
# Module 2 - Branch and Bound: one generic engine, three search orders
#
#   LIFO  (stack)     = depth-first B&B: dives fast, finds a first solution early
#   FIFO  (queue)     = breadth-first B&B: explores level by level
#   LC    (min-heap)  = least-cost / best-first B&B: always expands the node with
#                       the smallest lower bound => expands the fewest nodes
#
# The problem-specific part is supplied as closures:
#   root                  initial node (any R object, e.g. a list)
#   expand(node)          list of child nodes
#   bound(node)           LOWER bound on the cost of any completion (minimisation)
#   is_leaf(node)         TRUE for a complete solution
#   leaf_cost(node)       its cost
# A node is PRUNED when bound(node) >= best cost found so far.
# =============================================================================

#' Generic minimisation branch and bound.
#'
#' Problem: minimise cost over a tree of partial solutions.
#' Approach: keep a frontier of live nodes in a container whose discipline
#'   decides the strategy (stack / queue / heap keyed by bound). Pop a node,
#'   discard it if its bound cannot beat the incumbent, otherwise either record
#'   it (leaf) or push its children whose bound is still promising.
#' Time complexity: O(b^d) worst case (exponential), but the bound can prune
#'   most of the tree; LC expands the fewest nodes, LIFO needs the least memory.
#' Space complexity: O(frontier size): O(d*b) for LIFO, up to O(b^d) for FIFO/LC.
#'
#' @param strategy "lifo" | "fifo" | "lc"
#' @param ub initial upper bound (e.g. from a heuristic); Inf if none
#' @param max_nodes stop (and set `complete = FALSE`) after expanding this many
#' @return list(cost, node (best leaf), expanded, generated, pruned,
#'              max_frontier, complete)
branch_and_bound <- function(root, expand, bound, is_leaf, leaf_cost,
                             strategy = "lc", ub = Inf, max_nodes = Inf,
                             frontier_cap = 2000000L) {
  # node pool lives in an ENVIRONMENT: `pool[[i]] <- x` through <<- on a closure variable
  # would copy the whole list on every push (O(pool size) each -> quadratic).
  P <- new.env(parent = emptyenv()); P$pool <- vector("list", 1024L); P$n <- 0L
  new_node <- function(nd) {
    P$n <- P$n + 1L
    if (P$n > length(P$pool)) length(P$pool) <- 2L * length(P$pool)
    P$pool[[P$n]] <- nd
    P$n
  }
  st <- NULL; qu <- NULL; hp <- NULL
  if (strategy == "lifo") st <- stack_new(256L)
  else if (strategy == "fifo") qu <- queue_new(frontier_cap)
  else if (strategy == "lc") hp <- heap_new(1024L)
  else stop("unknown strategy: ", strategy)
  push <- function(id, key) {
    if (strategy == "lifo") stack_push(st, id)
    else if (strategy == "fifo") queue_push(qu, id)
    else heap_push(hp, key, id)
  }
  empty <- function() {
    if (strategy == "lifo") stack_empty(st) else if (strategy == "fifo") queue_empty(qu) else heap_empty(hp)
  }
  pop <- function() {
    if (strategy == "lifo") stack_pop(st) else if (strategy == "fifo") queue_pop(qu) else heap_pop(hp)$val
  }
  fsize <- function() {
    if (strategy == "lifo") stack_size(st) else if (strategy == "fifo") queue_size(qu) else heap_size(hp)
  }

  best <- ub; best_node <- NULL
  expanded <- 0L; generated <- 1L; pruned <- 0L; maxf <- 1L
  push(new_node(root), bound(root))
  complete <- TRUE
  while (!empty()) {
    id <- pop(); nd <- P$pool[[id]]; P$pool[id] <- list(NULL)   # free memory
    if (bound(nd) >= best) { pruned <- pruned + 1L; next }
    if (is_leaf(nd)) {
      c0 <- leaf_cost(nd)
      if (c0 < best) { best <- c0; best_node <- nd }
      next
    }
    if (expanded >= max_nodes) { complete <- FALSE; break }
    expanded <- expanded + 1L
    kids <- expand(nd)
    for (ch in kids) {
      generated <- generated + 1L
      b <- bound(ch)
      if (b < best) push(new_node(ch), b) else pruned <- pruned + 1L
    }
    if (fsize() > maxf) maxf <- fsize()
  }
  list(cost = best, node = best_node, expanded = expanded, generated = generated,
       pruned = pruned, max_frontier = maxf, complete = complete)
}

# ---------------------------------------------------------------------------
# Travelling salesman by branch and bound (reduced cost matrix bound)
# ---------------------------------------------------------------------------

#' Reduce a cost matrix: subtract each row's minimum, then each column's minimum
#' (only finite entries). Returns the reduced matrix and the total subtracted.
#' Time O(n^2).
tsp_reduce <- function(M) {
  n <- nrow(M); red <- 0
  for (i in seq_len(n)) {
    mn <- Inf; for (j in seq_len(n)) if (M[i, j] < mn) mn <- M[i, j]
    if (mn > 0 && mn < Inf) { for (j in seq_len(n)) if (M[i, j] < Inf) M[i, j] <- M[i, j] - mn; red <- red + mn }
  }
  for (j in seq_len(n)) {
    mn <- Inf; for (i in seq_len(n)) if (M[i, j] < mn) mn <- M[i, j]
    if (mn > 0 && mn < Inf) { for (i in seq_len(n)) if (M[i, j] < Inf) M[i, j] <- M[i, j] - mn; red <- red + mn }
  }
  list(M = M, red = red)
}

#' Travelling salesman by branch and bound.
#'
#' Problem: shortest closed tour visiting every city exactly once (asymmetric or
#'   symmetric costs D[i,j]); start city 1.
#' Approach: node = partial path from city 1. Lower bound = cost of the path so
#'   far + the "reduced-matrix" bound: every city must be left and entered
#'   exactly once, so the row/column minima of the remaining matrix must still
#'   be paid. Child (i -> j): copy the matrix, forbid row i and column j and the
#'   premature return j -> 1, re-reduce; bound = parent bound + D'[i,j] + the new
#'   reduction.
#' Time complexity: O(n!) worst case, O(n^2) work per node; the bound prunes
#'   most of the tree for metric instances.
#' Space complexity: O(n^2) per live node (matrix copy).
#' @return list(cost, tour (city order, starts at 1), expanded, generated, ...)
tsp_bb <- function(D, strategy = "lc", ub = Inf, max_nodes = Inf) {
  n <- nrow(D)
  if (n == 1L) return(list(cost = 0, tour = 1L, expanded = 0L, generated = 1L, complete = TRUE))
  M <- D; for (i in seq_len(n)) M[i, i] <- Inf
  r0 <- tsp_reduce(M)
  root <- list(M = r0$M, bound = r0$red, path = 1L, last = 1L)
  expand <- function(nd) {
    kids <- list(); i <- nd$last; len <- length(nd$path)
    visited <- logical(n); for (c in nd$path) visited[c] <- TRUE
    for (j in seq_len(n)) {
      if (visited[j]) next
      M2 <- nd$M
      step <- M2[i, j]
      for (k in seq_len(n)) { M2[i, k] <- Inf; M2[k, j] <- Inf }
      M2[j, 1L] <- Inf                       # forbid closing the tour early
      rr <- tsp_reduce(M2)
      kids[[length(kids) + 1L]] <- list(M = rr$M, bound = nd$bound + step + rr$red,
                                        path = c(nd$path, j), last = j)
    }
    kids
  }
  res <- branch_and_bound(root, expand, function(nd) nd$bound,
                          function(nd) length(nd$path) == n, function(nd) tour_cost(D, nd$path),
                          strategy = strategy, ub = ub, max_nodes = max_nodes)
  tour <- if (is.null(res$node)) integer(0) else res$node$path
  cost <- if (length(tour) == n) tour_cost(D, tour) else Inf
  list(cost = cost, tour = tour, expanded = res$expanded, generated = res$generated,
       pruned = res$pruned, max_frontier = res$max_frontier, complete = res$complete)
}

#' Length of the closed tour `tour` under cost matrix D. O(n).
tour_cost <- function(D, tour) {
  n <- length(tour); s <- 0
  for (i in seq_len(n)) s <- s + D[tour[i], tour[if (i == n) 1L else i + 1L]]
  s
}

#' REFERENCE: brute-force over all (n-1)! tours (Heap's algorithm, iterative).
#' Time O(n! * n), space O(n). n <= 10.
tsp_brute <- function(D) {
  n <- nrow(D)
  if (n <= 2L) return(list(cost = tour_cost(D, seq_len(n)), tour = seq_len(n)))
  perm <- 2:n; m <- n - 1L; cnt <- rep(1L, m)
  best <- tour_cost(D, c(1L, perm)); bt <- c(1L, perm)
  i <- 1L
  while (i <= m) {
    if (cnt[i] < i) {
      if (i %% 2L == 1L) { t <- perm[1L]; perm[1L] <- perm[i]; perm[i] <- t }
      else               { t <- perm[cnt[i]]; perm[cnt[i]] <- perm[i]; perm[i] <- t }
      cst <- tour_cost(D, c(1L, perm))
      if (cst < best) { best <- cst; bt <- c(1L, perm) }
      cnt[i] <- cnt[i] + 1L; i <- 1L
    } else { cnt[i] <- 1L; i <- i + 1L }
  }
  list(cost = best, tour = bt)
}
