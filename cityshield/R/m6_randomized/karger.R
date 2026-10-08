# =============================================================================
# Module 6 - Randomized algorithms: Karger's minimum cut
# CityShield use: find the WEAKEST cut of the road network - the smallest set of
#   road links whose loss disconnects the city (a resilience bottleneck).
# Instances are 50-150 node road subgraphs (see CLAUDE-level notes in report).
# =============================================================================

# ---- union-find (disjoint sets) ----------------------------------------------

#' Union-find with path compression and union by size. Nodes 1..n.
#' find / union: amortised O(alpha(n)) ~ O(1).  Space O(n).
uf_new <- function(n) {
  u <- new.env(parent = emptyenv())
  u$parent <- seq_len(n); u$size <- rep(1L, n); u$count <- n
  u
}
uf_find <- function(u, x) {
  r <- x
  while (u$parent[r] != r) r <- u$parent[r]
  while (u$parent[x] != r) { nx <- u$parent[x]; u$parent[x] <- r; x <- nx }   # path compression
  r
}
uf_union <- function(u, a, b) {
  ra <- uf_find(u, a); rb <- uf_find(u, b)
  if (ra == rb) return(FALSE)
  if (u$size[ra] < u$size[rb]) { t <- ra; ra <- rb; rb <- t }
  u$parent[rb] <- ra; u$size[ra] <- u$size[ra] + u$size[rb]; u$count <- u$count - 1L
  TRUE
}

# ---- Karger -----------------------------------------------------------------

#' One run of Karger's random-contraction algorithm.
#'
#' Problem: minimum cut of an undirected (multi)graph = fewest edges whose
#'   removal disconnects it.
#' Approach: repeat { pick a UNIFORMLY RANDOM remaining edge; contract its two
#'   endpoints into one super-node (union-find); drop self-loops } until only 2
#'   super-nodes remain. The edges between them are a cut.
#'   Self-loops are discarded lazily: a random edge whose endpoints are already
#'   merged is swapped out of the active prefix of the edge array.
#' Why it works: a fixed minimum cut C (size k) survives a contraction step
#'   unless an edge of C is picked; with n' super-nodes every node has degree >=
#'   k, so P(pick C-edge) <= 2/n'. Hence P(C survives all n-2 steps)
#'   >= prod (1 - 2/n') = 2 / (n (n-1)).
#' Time complexity per run: O(m * alpha(n)) amortised (each edge is removed at
#'   most once; each contraction is one union).
#' Space complexity: O(n + m).
#' @param n nodes; @param eu,ev edge endpoint vectors (multi-edges allowed)
#' @return list(cut = size of the cut found, side = logical per node,
#'              cut_edges = indices of the crossing edges)
karger_run <- function(n, eu, ev) {
  m <- length(eu)
  idx <- seq_len(m)                          # active edge ids occupy idx[1..active]
  uf <- uf_new(n)
  active <- m
  while (uf$count > 2L && active > 0L) {
    r <- sample.int(active, 1L)
    e <- idx[r]
    if (uf_find(uf, eu[e]) == uf_find(uf, ev[e])) {   # self-loop: remove it
      idx[r] <- idx[active]; idx[active] <- e; active <- active - 1L
      next
    }
    uf_union(uf, eu[e], ev[e])
  }
  # count crossing edges between the two remaining super-nodes
  roots <- integer(n); for (v in seq_len(n)) roots[v] <- uf_find(uf, v)
  first <- roots[1L]; side <- logical(n)
  for (v in seq_len(n)) side[v] <- roots[v] == first
  cut_edges <- integer(0)
  for (e in seq_len(m)) if (side[eu[e]] != side[ev[e]]) cut_edges <- c(cut_edges, e)
  list(cut = length(cut_edges), side = side, cut_edges = cut_edges)
}

#' Repeat Karger `trials` times and keep the best cut.
#' Time O(trials * m * alpha(n)).  If a min cut has probability >= p = 2/(n(n-1))
#' per trial, trials = ceil(ln(1/delta) / p) gives failure probability <= delta.
#' @return list(best = karger_run result, cuts = cut size of each trial)
karger_min_cut <- function(n, eu, ev, trials) {
  best <- NULL; cuts <- numeric(trials)
  for (t in seq_len(trials)) {
    r <- karger_run(n, eu, ev)
    cuts[t] <- r$cut
    if (is.null(best) || r$cut < best$cut) best <- r
  }
  list(best = best, cuts = cuts)
}

#' Trials needed so that Karger fails with probability <= delta (theory, using
#' the guaranteed per-trial success probability 2 / (n (n-1))). O(1).
karger_trials_needed <- function(n, delta = 0.01) ceiling(log(1 / delta) / (2 / (n * (n - 1))))

#' Guaranteed lower bound on the success probability after T trials:
#' 1 - (1 - 2/(n(n-1)))^T. O(1).
karger_success_bound <- function(n, trials) 1 - (1 - 2 / (n * (n - 1)))^trials

#' REFERENCE exact min cut by trying all 2^(n-1) bipartitions (n <= 16).
#' Time O(2^n * m).  Space O(n).
mincut_brute <- function(n, eu, ev) {
  best <- Inf
  for (mask in 1:(2^(n - 1) - 1)) {           # node n is always on side 0
    side <- logical(n); x <- mask
    for (v in seq_len(n - 1L)) { side[v] <- x %% 2 == 1; x <- x %/% 2 }
    c0 <- 0L
    for (e in seq_along(eu)) if (side[eu[e]] != side[ev[e]]) c0 <- c0 + 1L
    if (c0 < best) best <- c0
  }
  best
}

#' Undirected simple edge list (one entry per road LINK) from a CSR graph:
#' u->v and v->u collapse to one edge, self loops and duplicates dropped.
#' O(V + E) using a mark array per node.
csr_undirected_edges <- function(g) {
  adj <- csr_undirected_adj(g)
  eu <- integer(0); ev <- integer(0)
  for (u in seq_len(g$n)) for (v in adj[[u]]) if (u < v) { eu <- c(eu, u); ev <- c(ev, v) }
  list(n = g$n, eu = eu, ev = ev)
}

#' Is the undirected graph connected? Union-find over the edges. O(E alpha(V)).
edges_connected <- function(n, eu, ev) {
  if (n <= 1L) return(TRUE)
  uf <- uf_new(n); for (e in seq_along(eu)) uf_union(uf, eu[e], ev[e])
  uf$count == 1L
}
