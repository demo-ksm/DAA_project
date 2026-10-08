# =============================================================================
# CityShield application layer for Module 4 (thin wrappers, no new algorithms)
#   - evacuation_network():  road graph + supplies + shelters -> flow network
#   - dispatch_feasibility(): which ambulance can reach which incident in time
# =============================================================================

#' Build the evacuation flow network on (a subgraph of) the road network.
#'
#' Super-source S -> each danger node (capacity = people to evacuate),
#' road edges keep their capacity (vehicles/h * persons per vehicle),
#' each shelter node -> super-sink T (capacity = shelter space).
#' CSR edge k corresponds to flow edge slot 2k-1 (see flow_from_csr).
#' Time O(V + E).
#' @return list(fg, s, t, n_road_edges, danger, shelters)
evacuation_network <- function(g, danger, people, shelters, shelter_cap,
                               persons_per_vehicle = 3) {
  cap <- g$extra$cap * persons_per_vehicle
  fg <- flow_new(g$n + 2L, g$m + length(danger) + length(shelters))
  for (u in seq_len(g$n)) {
    lo <- g$offset[u]
    for (k in seq_len(g$offset[u + 1L] - lo)) {
      e <- lo + k - 1L
      flow_add_edge(fg, u, g$target[e], cap[e])
    }
  }
  s <- g$n + 1L; t <- g$n + 2L
  for (i in seq_along(danger))   flow_add_edge(fg, s, danger[i], people[i])
  for (j in seq_along(shelters)) flow_add_edge(fg, shelters[j], t, shelter_cap[j])
  list(fg = fg, s = s, t = t, n_road_edges = g$m, danger = danger, shelters = shelters)
}

#' Flow on every CSR road edge after a max-flow run on an evacuation network.
#' O(E).
evacuation_edge_flow <- function(net) {
  out <- numeric(net$n_road_edges)
  for (k in seq_len(net$n_road_edges)) {
    e <- 2L * k - 1L
    out[k] <- net$fg$cap0[e] - net$fg$cap[e]
  }
  out
}

#' Compatibility pairs (ambulance i, incident j) where the travel time from the
#' ambulance's station is within `max_time` seconds.
#' One Dijkstra per ambulance, stopped at `max_time` (limit).
#' Time: O(A * (V' + E') log V') where V', E' = nodes/edges within the radius.
#' @param amb_nodes node of each ambulance; @param inc_nodes node of each incident
#' @return list(l, r, time) edge lists for max_bipartite_matching
dispatch_feasibility <- function(g, amb_nodes, inc_nodes, max_time) {
  l <- integer(0); r <- integer(0); tt <- numeric(0)
  for (i in seq_along(amb_nodes)) {
    d <- dijkstra(g, amb_nodes[i], limit = max_time)$dist
    for (j in seq_along(inc_nodes)) {
      if (d[inc_nodes[j]] <= max_time) {
        l <- c(l, i); r <- c(r, j); tt <- c(tt, d[inc_nodes[j]])
      }
    }
  }
  list(l = l, r = r, time = tt)
}
