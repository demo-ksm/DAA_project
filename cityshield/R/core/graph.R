# =============================================================================
# core/graph.R  -  Custom graph representations for CityShield
#
# Two representations (design choice approved in Phase 0):
#
#   1. CSR  (Compressed Sparse Row)  - for shortest-path style algorithms.
#        offset : integer[V+1]  neighbours of u are target[offset[u] .. offset[u+1]-1]
#        target : integer[E]    head node of each directed edge
#        weight : numeric[E]    edge weight (travel time in seconds)
#        x, y   : numeric[V]    node coordinates (lon, lat)
#      Neighbour scan of u costs O(deg(u)); total memory O(V + E).
#
#   2. Residual graph - for flow algorithms.
#        Every edge is stored with its reverse edge as the pair (2k-1, 2k); see
#        rev_edge(). Adjacency uses linked lists (head/nxt) so edges append in
#        O(1). Held in an environment because R copies vectors on modify.
#
# Nothing here uses igraph/sf: those are only used (elsewhere) to LOAD the OSM
# data and to VALIDATE results in tests.
# =============================================================================

# ---- 1. CSR graph -----------------------------------------------------------

#' Build a CSR graph from an edge list (counting sort by source node).
#'
#' @description
#' Problem: convert an unsorted directed edge list into a structure that lets us
#' enumerate the out-neighbours of any node in O(deg).
#'
#' Approach: counting sort on the source id.
#'   pass 1: count out-degree of every node
#'   pass 2: prefix sums give `offset`
#'   pass 3: drop every edge into its slot
#'
#' Time complexity: O(V + E).  Space complexity: O(V + E).
#'
#' @param n      number of nodes (ids are 1..n)
#' @param from,to  integer vectors of length E (edge endpoints)
#' @param w      numeric vector of length E (edge weights)
#' @param x,y    node coordinates (length n), may be NULL
#' @param extra  optional named list of per-edge vectors (input order), e.g.
#'   list(cap = ..., hw = ...); stored permuted in `g$extra`
#' @return list of class "csr_graph"; `eid` lets callers permute any extra
#'   per-edge attribute: attr_in_csr_order <- attr_in_input_order[g$eid]
csr_build <- function(n, from, to, w, x = NULL, y = NULL, extra = NULL) {
  m <- length(from)
  stopifnot(length(to) == m, length(w) == m)

  # pass 1: out-degree counts
  deg <- integer(n)
  for (e in seq_len(m)) deg[from[e]] <- deg[from[e]] + 1L

  # pass 2: prefix sums -> offset (1-based, offset[n+1] = m + 1)
  offset <- integer(n + 1L)
  offset[1L] <- 1L
  for (u in seq_len(n)) offset[u + 1L] <- offset[u] + deg[u]

  # pass 3: place each edge in its slot ("cursor" = next free slot of u)
  cursor <- offset[seq_len(n)]
  target <- integer(m)
  weight <- numeric(m)
  eid    <- integer(m)          # eid[slot] = index of that edge in the INPUT lists
  for (e in seq_len(m)) {
    u <- from[e]
    slot <- cursor[u]
    target[slot] <- to[e]
    weight[slot] <- w[e]
    eid[slot] <- e
    cursor[u] <- slot + 1L
  }

  # per-edge extras (capacity, road class, ...) are permuted to CSR order
  ex <- list()
  for (nm in names(extra)) ex[[nm]] <- extra[[nm]][eid]

  structure(
    list(n = n, m = m, offset = offset, target = target, weight = weight,
         eid = eid, x = x, y = y, extra = ex),
    class = "csr_graph"
  )
}

#' Out-degree of node u in a CSR graph. O(1).
csr_degree <- function(g, u) g$offset[u + 1L] - g$offset[u]

#' Reverse a CSR graph (swap edge direction). O(V + E).
#' Needed e.g. for "distance TO a hospital" queries (run Dijkstra on the reverse).
csr_reverse <- function(g) {
  from <- integer(g$m)
  for (u in seq_len(g$n)) {
    lo <- g$offset[u]
    hi <- g$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) from[e] <- u
  }
  # extras are already in CSR order of g, and edge k of the new input list is
  # CSR edge k of g, so they pass through unchanged.
  csr_build(g$n, from = g$target, to = from, w = g$weight, x = g$x, y = g$y,
            extra = g$extra)
}

#' Extract the node-induced subgraph on the first `k` nodes of a BFS from `root`.
#'
#' Used to produce the 50-150 node graphs required by TSP B&B, Karger and exact
#' covers, and the growing sizes used in benchmarks.
#'
#' Time: O(V + E).  Space: O(V).
#' @return csr_graph with nodes renumbered 1..k and attr "orig_id"
csr_bfs_subgraph <- function(g, root, k) {
  k <- min(k, g$n)
  new_id <- integer(g$n)             # 0 = not selected
  orig <- integer(k)
  q <- queue_new(g$n)                # defined in core/structures.R
  queue_push(q, root)
  new_id[root] <- 1L
  orig[1L] <- root
  count <- 1L
  while (!queue_empty(q) && count < k) {
    u <- queue_pop(q)
    lo <- g$offset[u]; hi <- g$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) {
      v <- g$target[e]
      if (new_id[v] == 0L && count < k) {
        count <- count + 1L
        new_id[v] <- count
        orig[count] <- v
        queue_push(q, v)
      }
    }
  }
  # collect edges whose both endpoints are selected
  cap <- 1024L; fr <- integer(cap); to <- integer(cap); w <- numeric(cap)
  src <- integer(cap)               # CSR edge id in g, to carry extras over
  cnt <- 0L
  for (i in seq_len(count)) {
    u <- orig[i]
    lo <- g$offset[u]; hi <- g$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) {
      v <- new_id[g$target[e]]
      if (v != 0L) {
        cnt <- cnt + 1L
        if (cnt > cap) {             # amortised doubling: O(1) per append
          cap <- cap * 2L
          length(fr) <- cap; length(to) <- cap; length(w) <- cap; length(src) <- cap
        }
        fr[cnt] <- i; to[cnt] <- v; w[cnt] <- g$weight[e]; src[cnt] <- e
      }
    }
  }
  ex <- list()
  for (nm in names(g$extra)) ex[[nm]] <- g$extra[[nm]][src[seq_len(cnt)]]
  sub <- csr_build(count, fr[seq_len(cnt)], to[seq_len(cnt)], w[seq_len(cnt)],
                   extra = ex,
                   x = if (!is.null(g$x)) g$x[orig[seq_len(count)]] else NULL,
                   y = if (!is.null(g$y)) g$y[orig[seq_len(count)]] else NULL)
  attr(sub, "orig_id") <- orig[seq_len(count)]
  sub
}

#' Undirected adjacency lists from a CSR graph (edges treated as undirected,
#' self loops dropped, duplicates removed). Time O(V + E) with a mark array.
csr_undirected_adj <- function(g) {
  n <- g$n
  nb <- vector("list", n)
  for (u in seq_len(n)) nb[[u]] <- integer(0)
  mark <- integer(n)
  for (u in seq_len(n)) {
    lo <- g$offset[u]
    for (k in seq_len(g$offset[u + 1L] - lo)) {
      v <- g$target[lo + k - 1L]
      if (v != u) {
        nb[[u]] <- c(nb[[u]], v); nb[[v]] <- c(nb[[v]], u)
      }
    }
  }
  for (u in seq_len(n)) {                    # remove duplicates
    x <- nb[[u]]; keep <- integer(length(x)); m <- 0L
    for (v in x) if (mark[v] != u) { mark[v] <- u; m <- m + 1L; keep[m] <- v }
    nb[[u]] <- keep[seq_len(m)]
  }
  nb
}

# ---- 2. Residual graph for flow --------------------------------------------
#
# R copies vectors when they are modified inside a function ("copy-on-modify"),
# so a list-based graph would cost O(E) per added edge. We therefore keep the
# flow graph in an ENVIRONMENT (reference semantics): updates are in place, O(1).

#' Create an empty residual graph with room for `max_edges` logical edges
#' (each uses TWO slots: forward + reverse, stored as the pair (2k-1, 2k)).
#'
#' Adjacency lists are linked lists (`head` / `nxt`) so adding an edge is O(1).
#' @return environment of class "flow_graph"
flow_new <- function(n, max_edges) {
  cap2 <- 2L * as.integer(max_edges)
  fg <- new.env(parent = emptyenv())
  fg$n    <- n
  fg$head <- integer(n)        # head[u] = first edge id out of u (0 = none)
  fg$nxt  <- integer(cap2)     # next edge in u's list
  fg$to   <- integer(cap2)     # endpoint of each edge slot
  fg$cap  <- numeric(cap2)     # residual capacity (mutated during a run)
  fg$cap0 <- numeric(cap2)     # original capacity (to reset / read the flow)
  fg$ne   <- 0L                # number of edge slots used
  class(fg) <- "flow_graph"
  fg
}

#' Partner (reverse) edge index. Edges come in pairs (1,2), (3,4), ...
rev_edge <- function(e) if (e %% 2L == 1L) e + 1L else e - 1L

#' Add edge u->v with capacity c (plus reverse edge v->u with capacity 0).
#' Modifies `fg` in place. Time O(1).  Returns the id of the forward edge.
flow_add_edge <- function(fg, u, v, c) {
  e <- fg$ne + 1L
  fg$to[e] <- v;  fg$cap[e] <- c;  fg$cap0[e] <- c
  fg$nxt[e] <- fg$head[u]; fg$head[u] <- e
  e2 <- e + 1L
  fg$to[e2] <- u; fg$cap[e2] <- 0; fg$cap0[e2] <- 0
  fg$nxt[e2] <- fg$head[v]; fg$head[v] <- e2
  fg$ne <- e2
  e
}

#' Restore all residual capacities to the original ones. O(E).
flow_reset <- function(fg) {
  for (e in seq_len(fg$ne)) fg$cap[e] <- fg$cap0[e]
  invisible(fg)
}

#' Build a flow graph from a CSR graph; `capacity[e]` is the capacity of CSR
#' edge e (e.g. road capacity derived from road class). Time O(V + E).
flow_from_csr <- function(g, capacity) {
  fg <- flow_new(g$n, g$m)
  for (u in seq_len(g$n)) {
    lo <- g$offset[u]; hi <- g$offset[u + 1L] - 1L
    if (hi >= lo) for (e in lo:hi) flow_add_edge(fg, u, g$target[e], capacity[e])
  }
  fg
}
