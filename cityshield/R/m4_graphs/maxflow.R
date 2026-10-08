# =============================================================================
# Module 4 - Maximum flow: Push-Relabel
# CityShield use: evacuation planning (people/vehicles per hour from danger
#                 zones to shelters through roads with limited capacity) and
#                 ambulance-to-incident matching (via bipartite_matching.R).
#
# Works on a `flow_graph` (core/graph.R): edge pairs (e, rev_edge(e)),
# residual capacities in fg$cap, originals in fg$cap0. It MODIFIES fg in place;
# call flow_reset(fg) to run again on the same network.
# Flow on forward edge e after a run = fg$cap0[e] - fg$cap[e].
# =============================================================================

#' Exact distance-to-sink labels by reverse BFS on the residual graph.
#' Nodes that cannot reach t get label n. O(V + E).
.pr_global_relabel <- function(fg, s, t, height) {
  n <- fg$n
  for (i in seq_len(n)) height[i] <- n
  height[t] <- 0L
  q <- queue_new(n)
  queue_push(q, t)
  while (!queue_empty(q)) {
    v <- queue_pop(q)
    e <- fg$head[v]
    while (e != 0L) {
      u <- fg$to[e]                         # edge e: v -> u, so rev_edge(e): u -> v
      if (u != s && height[u] == n && fg$cap[rev_edge(e)] > 0) {
        height[u] <- height[v] + 1L
        queue_push(q, u)
      }
      e <- fg$nxt[e]
    }
  }
  height[s] <- n
  height
}

#' Push-Relabel (FIFO selection) max flow.
#'
#' Problem: maximum s-t flow.
#' Approach: instead of whole augmenting paths, keep a PREFLOW (nodes may hold
#'   excess) and local heights. Saturate all source edges, then repeatedly
#'   take an active node (excess > 0) from a FIFO queue and `discharge` it:
#'   push excess along admissible edges (residual > 0 and height[u] ==
#'   height[v] + 1); if stuck, RELABEL to 1 + the lowest residual neighbour.
#'   A current-edge pointer avoids rescanning edges. Initial exact labels come
#'   from a reverse BFS (global relabeling); the GAP heuristic lifts every node
#'   above an emptied height level straight to n+1 (none of them can reach t).
#'   PHASE 1 (default, `complete = FALSE`): a node whose height reaches n can no
#'   longer reach the sink, so its excess is left alone - when no other active
#'   node remains, excess[t] IS the maximum flow value (the stranded excess only
#'   needs to travel back to s to turn the preflow into a flow). Phase 2
#'   (`complete = TRUE`) returns it, giving a valid edge flow (needed for
#'   matchings and flow decomposition).
#' Time complexity: O(V^3) for FIFO selection (O(V^2 E) for generic).
#' Space complexity: O(V) extra.
#'
#' @param complete TRUE = also return stranded excess to the source (valid flow)
#' @param gap use the gap heuristic
#' @return list(value, pushes, relabels)
push_relabel <- function(fg, s, t, global_relabel = TRUE, complete = FALSE, gap = TRUE) {
  n <- fg$n
  excess <- numeric(n)
  height <- integer(n)
  cur <- fg$head                             # current-arc pointers
  in_q <- logical(n)
  q <- queue_new(n)
  pushes <- 0L; relabels <- 0L

  if (global_relabel) height <- .pr_global_relabel(fg, s, t, height) else height[s] <- n
  cnt <- integer(2L * n + 2L)                # cnt[h+1] = number of non-source nodes at height h
  for (v in seq_len(n)) if (v != s) cnt[height[v] + 1L] <- cnt[height[v] + 1L] + 1L
  # saturate every edge out of the source
  e <- fg$head[s]
  while (e != 0L) {
    c <- fg$cap[e]
    if (c > 0) {
      v <- fg$to[e]
      fg$cap[e] <- 0; r <- rev_edge(e); fg$cap[r] <- fg$cap[r] + c
      excess[v] <- excess[v] + c; excess[s] <- excess[s] - c
      if (v != s && v != t && !in_q[v]) { queue_push(q, v); in_q[v] <- TRUE }
    }
    e <- fg$nxt[e]
  }

  while (!queue_empty(q)) {
    u <- queue_pop(q); in_q[u] <- FALSE
    if (!complete && height[u] >= n) next    # phase 1: u can no longer reach t
    # discharge u
    while (excess[u] > 0) {
      e <- cur[u]
      if (e == 0L) {                         # no admissible edge left: relabel
        mh <- 2L * n
        e2 <- fg$head[u]
        while (e2 != 0L) {
          if (fg$cap[e2] > 0 && height[fg$to[e2]] < mh) mh <- height[fg$to[e2]]
          e2 <- fg$nxt[e2]
        }
        old <- height[u]
        cnt[old + 1L] <- cnt[old + 1L] - 1L
        height[u] <- mh + 1L
        cnt[height[u] + 1L] <- cnt[height[u] + 1L] + 1L
        relabels <- relabels + 1L
        cur[u] <- fg$head[u]
        if (gap && old < n && cnt[old + 1L] == 0L) {         # gap: levels above `old` are cut off from t
          for (v in seq_len(n)) {
            if (v != s && height[v] > old && height[v] < n + 1L) {
              cnt[height[v] + 1L] <- cnt[height[v] + 1L] - 1L
              height[v] <- n + 1L; cnt[n + 2L] <- cnt[n + 2L] + 1L
              cur[v] <- fg$head[v]
            }
          }
        }
        if (!complete && height[u] >= n) break               # phase 1: stop working on u
        next
      }
      v <- fg$to[e]
      if (fg$cap[e] > 0 && height[u] == height[v] + 1L) {
        d <- if (excess[u] < fg$cap[e]) excess[u] else fg$cap[e]
        fg$cap[e] <- fg$cap[e] - d
        r <- rev_edge(e); fg$cap[r] <- fg$cap[r] + d
        excess[u] <- excess[u] - d; excess[v] <- excess[v] + d
        pushes <- pushes + 1L
        if (v != s && v != t && !in_q[v]) { queue_push(q, v); in_q[v] <- TRUE }
      } else {
        cur[u] <- fg$nxt[e]
      }
    }
  }
  list(value = excess[t], pushes = pushes, relabels = relabels)
}

#' Flow carried by each forward edge slot after a run: cap0 - cap (>= 0 only
#' for forward edges; reverse slots return their negative). O(E).
flow_values <- function(fg) {
  out <- numeric(fg$ne)
  for (e in seq_len(fg$ne)) out[e] <- fg$cap0[e] - fg$cap[e]
  out
}

#' Minimum s-t cut from a finished max-flow run: nodes reachable from s in the
#' residual graph form the source side. Returns logical vector (TRUE = source
#' side). O(V + E). Used to verify max-flow = min-cut in tests.
min_cut_side <- function(fg, s) {
  seen <- logical(fg$n)
  q <- queue_new(fg$n)
  queue_push(q, s); seen[s] <- TRUE
  while (!queue_empty(q)) {
    u <- queue_pop(q)
    e <- fg$head[u]
    while (e != 0L) {
      v <- fg$to[e]
      if (fg$cap[e] > 0 && !seen[v]) { seen[v] <- TRUE; queue_push(q, v) }
      e <- fg$nxt[e]
    }
  }
  seen
}
