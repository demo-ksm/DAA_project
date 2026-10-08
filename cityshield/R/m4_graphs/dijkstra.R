# =============================================================================
# Module 4 (supporting) - Dijkstra's algorithm
# Used by: dispatch (nearest station), ambulance feasibility edges, TSP
# distance matrices, coverage/set-cover instances.
# =============================================================================

#' Dijkstra single/multi-source shortest paths (binary-heap version).
#'
#' Problem: shortest travel time from one node (or the nearest of several
#'   source nodes) to every node of a graph with NON-NEGATIVE weights.
#' Approach: greedy. Repeatedly settle the unsettled node with the smallest
#'   tentative distance (popped from a min-heap) and relax its out-edges.
#'   Stale heap entries are skipped ("lazy deletion"). Passing several
#'   `sources` starts them all at distance 0 = nearest-facility distances.
#' Time complexity: O((V + E) log V) with a binary heap.
#' Space complexity: O(V + E) (dist, parent, heap).
#'
#' @param g       csr_graph
#' @param sources integer vector of start nodes
#' @param target  optional node; stop as soon as it is settled
#' @param limit   stop expanding nodes farther than this distance
#' @return list(dist, parent, origin = index (in `sources`) of the source each node
#'   was reached from = "nearest facility" label, settled = number of settled nodes)
dijkstra <- function(g, sources, target = 0L, limit = Inf) {
  n <- g$n
  dist <- rep(Inf, n)
  parent <- integer(n)
  origin <- integer(n)
  done <- logical(n)
  h <- heap_new(1024L)
  for (i in seq_along(sources)) { s <- sources[i]; dist[s] <- 0; origin[s] <- i; heap_push(h, 0, s) }
  settled <- 0L
  while (!heap_empty(h)) {
    top <- heap_pop(h)
    u <- top$val
    if (done[u]) next                       # stale entry
    if (top$key > limit) break
    done[u] <- TRUE
    settled <- settled + 1L
    if (u == target) break
    lo <- g$offset[u]
    for (k in seq_len(g$offset[u + 1L] - lo)) {
      e <- lo + k - 1L
      v <- g$target[e]
      nd <- dist[u] + g$weight[e]
      if (nd < dist[v]) {
        dist[v] <- nd
        parent[v] <- u
        origin[v] <- origin[u]
        heap_push(h, nd, v)
      }
    }
  }
  list(dist = dist, parent = parent, origin = origin, settled = settled)
}

#' Dijkstra with an unsorted array instead of a heap (REFERENCE implementation).
#'
#' Problem: same as dijkstra(); used as the brute-force baseline.
#' Approach: pick the minimum unsettled node by a linear scan.
#' Time complexity: O(V^2 + E).  Space complexity: O(V).
dijkstra_array <- function(g, source) {
  n <- g$n
  dist <- rep(Inf, n); parent <- integer(n); done <- logical(n)
  dist[source] <- 0
  for (iter in seq_len(n)) {
    u <- 0L; best <- Inf
    for (i in seq_len(n)) if (!done[i] && dist[i] < best) { best <- dist[i]; u <- i }
    if (u == 0L) break
    done[u] <- TRUE
    lo <- g$offset[u]
    for (k in seq_len(g$offset[u + 1L] - lo)) {
      e <- lo + k - 1L; v <- g$target[e]
      if (dist[u] + g$weight[e] < dist[v]) { dist[v] <- dist[u] + g$weight[e]; parent[v] <- u }
    }
  }
  list(dist = dist, parent = parent)
}

#' Reconstruct the node path src -> dst from a parent vector. O(path length).
#' Returns integer(0) when dst is unreachable.
path_from_parent <- function(parent, src, dst, dist = NULL) {
  if (!is.null(dist) && dist[dst] == Inf) return(integer(0))
  len <- 1L; v <- dst
  while (v != src) {
    v <- parent[v]
    if (v == 0L) return(integer(0))
    len <- len + 1L
  }
  path <- integer(len); v <- dst
  for (i in len:1) { path[i] <- v; v <- parent[v] }
  path
}
