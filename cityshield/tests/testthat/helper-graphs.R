# Shared helpers for tests (testthat sources helper*.R automatically).

# Random directed multigraph as plain edge lists.
rand_edges <- function(n, m, seed, wmin = 1, wmax = 10, integer_w = FALSE) {
  set.seed(seed)
  w <- runif(m, wmin, wmax)
  if (integer_w) w <- round(w)
  list(n = n, from = sample.int(n, m, replace = TRUE),
       to = sample.int(n, m, replace = TRUE), w = w)
}

rand_csr <- function(n, m, seed, ...) {
  el <- rand_edges(n, m, seed, ...)
  csr_build(el$n, el$from, el$to, el$w)
}

# igraph copy of an edge list (for VALIDATION only).
to_igraph <- function(el, weights = TRUE) {
  ig <- igraph::make_graph(rbind(el$from, el$to), n = el$n, directed = TRUE)
  if (weights) igraph::E(ig)$weight <- el$w
  ig
}

# CSR graph back to an igraph (parallel edges kept), weights = travel time.
csr_to_igraph <- function(g) {
  from <- integer(g$m)
  for (u in seq_len(g$n)) {
    lo <- g$offset[u]
    for (k in seq_len(g$offset[u + 1L] - lo)) from[lo + k - 1L] <- u
  }
  ig <- igraph::make_graph(rbind(from, g$target), n = g$n, directed = TRUE)
  igraph::E(ig)$weight <- g$weight
  ig
}

# Flow network (edge list with integer capacities) for max-flow tests.
rand_flow_edges <- function(n, m, seed, cmax = 20) {
  set.seed(seed)
  list(n = n, from = sample.int(n, m, replace = TRUE),
       to = sample.int(n, m, replace = TRUE),
       cap = sample.int(cmax, m, replace = TRUE))
}

build_fg <- function(fe) {
  fg <- flow_new(fe$n, length(fe$from))
  for (i in seq_along(fe$from)) flow_add_edge(fg, fe$from[i], fe$to[i], fe$cap[i])
  fg
}
