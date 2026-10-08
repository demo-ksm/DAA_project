# =============================================================================
# core/study_area.R - helpers that cut the city graph down to the sizes each
# algorithm can handle in R (used by the precompute script, the dashboard and
# the benchmarks). Plain loops; none of this is a graded algorithm.
# =============================================================================

#' Metres east/north of a reference point (equirectangular approximation, fine
#' at city scale). Vectorised on purpose: this is coordinate preprocessing.
local_xy <- function(lon, lat, lon0, lat0) {
  list(x = (lon - lon0) * 111320 * cos(lat0 * pi / 180), y = (lat - lat0) * 110574)
}

#' Node-induced subgraph for the nodes where `keep` is TRUE (renumbered 1..k,
#' attr "orig_id" maps back). Keeps coordinates and per-edge extras. O(V + E).
csr_induced <- function(g, keep) {
  new_id <- integer(g$n); k <- 0L
  for (u in seq_len(g$n)) if (keep[u]) { k <- k + 1L; new_id[u] <- k }
  orig <- integer(k); for (u in seq_len(g$n)) if (keep[u]) orig[new_id[u]] <- u
  cap <- 1024L; fr <- integer(cap); to <- integer(cap); w <- numeric(cap); src <- integer(cap); cnt <- 0L
  for (u in seq_len(g$n)) {
    if (!keep[u]) next
    lo <- g$offset[u]
    for (j in seq_len(g$offset[u + 1L] - lo)) {
      e <- lo + j - 1L; v <- new_id[g$target[e]]
      if (v != 0L) {
        cnt <- cnt + 1L
        if (cnt > cap) { cap <- 2L * cap; length(fr) <- cap; length(to) <- cap; length(w) <- cap; length(src) <- cap }
        fr[cnt] <- new_id[u]; to[cnt] <- v; w[cnt] <- g$weight[e]; src[cnt] <- e
      }
    }
  }
  ex <- list(); for (nm in names(g$extra)) ex[[nm]] <- g$extra[[nm]][src[seq_len(cnt)]]
  sub <- csr_build(k, fr[seq_len(cnt)], to[seq_len(cnt)], w[seq_len(cnt)], extra = ex,
                   x = g$x[orig], y = g$y[orig])
  attr(sub, "orig_id") <- orig
  sub
}

#' Symmetrised shortest-path metric of a small graph: one Dijkstra run from every
#' node, then (D + t(D))/2 - still a metric, and symmetric as the TSP
#' 2-approximation needs. O(V (V + E) log V).
metric_closure <- function(g) {
  n <- g$n
  D <- matrix(Inf, n, n)
  for (s in seq_len(n)) D[s, ] <- dijkstra(g, s)$dist
  S <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) S[i, j] <- (D[i, j] + D[j, i]) / 2
  S
}

#' A k-node BFS subgraph of g whose pairwise distances are all finite (strongly
#' connected enough for metric algorithms). Tries up to `tries` random roots.
#' @return list(g, D)
study_instance <- function(g, k, seed, tries = 60L) {
  for (s in seq_len(tries)) {
    set.seed(seed * 1000L + s)
    sub <- csr_bfs_subgraph(g, sample.int(g$n, 1L), k)
    if (sub$n < k) next
    D <- metric_closure(sub)
    ok <- TRUE
    for (i in seq_len(sub$n)) for (j in seq_len(sub$n)) if (!is.finite(D[i, j])) { ok <- FALSE; break }
    if (ok) return(list(g = sub, D = D))
  }
  stop("no strongly connected ", k, "-node subgraph found")
}

#' Nodes within `radius_m` metres of (lon0, lat0) (logical vector). O(V).
nodes_within <- function(g, lon0, lat0, radius_m) {
  p <- local_xy(g$x, g$y, lon0, lat0)
  keep <- logical(g$n)
  for (i in seq_len(g$n)) keep[i] <- p$x[i]^2 + p$y[i]^2 <= radius_m^2
  keep
}
