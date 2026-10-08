# Phase 6 - Precompute everything expensive for the dashboard.
#
#   Rscript scripts/05_precompute.R        (CITYSHIELD_SYNTH=1 for the dev grid)
#
# Writes data/precomputed/*.rds. The Shiny app only READS these files, so every
# tab opens instantly. Cheap algorithms (knapsacks, small TSP, Karger, KMP search,
# quicksort ranking, matching ...) run LIVE in the app instead.
#
# Sizes follow the project rules: the whole city graph for Dijkstra/Dispatch and
# the local evacuation windows; 100-150 node study subgraphs for the NP-hard
# placement, Karger and TSP instances.

source("R/core/load_all.R")
g <- load_chennai_graph()
syn <- load_synthetic()
inc <- syn$incidents; st <- syn$stations; hos <- syn$hospitals
out_dir <- data_path("precomputed"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
save_pc <- function(obj, name) { saveRDS(obj, file.path(out_dir, paste0(name, ".rds"))); cat("saved", name, "\n"); flush.console() }
cat(sprintf("graph: %d nodes, %d edges\n", g$n, g$m))
set.seed(2026)

# =============================================================================
# 1. DISPATCH  (Dijkstra, quicksort ranking; the matching runs live in the app)
# =============================================================================
dj_st <- dijkstra(g, st$node)                       # distance FROM the nearest station
rg <- csr_reverse(g)
dj_h <- dijkstra(rg, hos$node)                       # distance FROM a node TO its nearest hospital
inc$station_time <- dj_st$dist[inc$node]; inc$station <- st$id[dj_st$origin[inc$node]]
inc$hospital_time <- dj_h$dist[inc$node]; inc$hospital <- hos$id[dj_h$origin[inc$node]]

# surge-hour batch: the 25 most urgent incidents of the storm window, matched to ambulances
storm <- which(inc$time_min >= 8 * 1440 & inc$time_min < 10 * 1440)
ord <- rank_incidents(inc$severity[storm], inc$demand[storm])
batch <- storm[ord[seq_len(min(25L, length(ord)))]]
amb <- st[st$type == "ambulance", ]
times <- matrix(Inf, nrow(amb), length(batch))
for (a in seq_len(nrow(amb))) {
  d <- dijkstra(g, amb$node[a], limit = 3600)$dist
  for (j in seq_along(batch)) times[a, j] <- d[inc$node[batch[j]]]
}
save_pc(list(incidents = inc, batch = batch, amb = amb, times = times), "dispatch")

# =============================================================================
# 2. EVACUATION  (push-relabel max flow, sweep line, Graham scan)
# =============================================================================
# pick 3 well separated hotspots by incident density
cand <- sample.int(nrow(inc), min(200L, nrow(inc)))
dens <- integer(length(cand))
for (i in seq_along(cand)) {
  p <- local_xy(inc$lon, inc$lat, inc$lon[cand[i]], inc$lat[cand[i]])
  dens[i] <- sum(p$x^2 + p$y^2 <= 1500^2)
}
centres <- integer(0)
for (i in merge_sort_perm(-dens)) {
  ok <- TRUE
  for (c0 in centres) { p <- local_xy(inc$lon[cand[i]], inc$lat[cand[i]], inc$lon[c0], inc$lat[c0]); if (sqrt(p$x^2 + p$y^2) < 4000) ok <- FALSE }
  if (ok) centres <- c(centres, cand[i])
  if (length(centres) == 3L) break
}
R_FLOOD <- 1000; R_WIN <- 2600
scen <- list()
for (s in seq_along(centres)) {
  lon0 <- inc$lon[centres[s]]; lat0 <- inc$lat[centres[s]]
  sub <- csr_induced(g, nodes_within(g, lon0, lat0, R_WIN))
  p <- local_xy(sub$x, sub$y, lon0, lat0); d <- sqrt(p$x^2 + p$y^2)
  dang <- integer(0); for (i in seq_len(sub$n)) if (d[i] <= R_FLOOD) dang <- c(dang, i)
  ring <- integer(0); for (i in seq_len(sub$n)) if (d[i] >= 1900 && d[i] <= 2500) ring <- c(ring, i)
  shel <- ring[sample.int(length(ring), min(10L, length(ring)))]
  people <- rep(25, length(dang)); scap <- rep(ceiling(sum(people) / 5), length(shel))
  net <- evacuation_network(sub, dang, people, shel, scap)
  pr <- push_relabel(net$fg, net$s, net$t, complete = TRUE)   # complete: a valid flow, so edge flows can be drawn
  flows <- evacuation_edge_flow(net)
  # --- geometry: which roads cross the flood boundary? (sweep line)
  ne <- 0L; er <- integer(0)
  for (u in seq_len(sub$n)) for (k in seq_len(sub$offset[u + 1L] - sub$offset[u])) {
    e <- sub$offset[u] + k - 1L; v <- sub$target[e]
    has_rev <- FALSE
    if (v < u) for (k2 in seq_len(sub$offset[v + 1L] - sub$offset[v])) if (sub$target[sub$offset[v] + k2 - 1L] == u) { has_rev <- TRUE; break }
    if (u < v || !has_rev) { ne <- ne + 1L; er <- c(er, e) }
  }
  from_of <- integer(sub$m); for (u in seq_len(sub$n)) for (k in seq_len(sub$offset[u + 1L] - sub$offset[u])) from_of[sub$offset[u] + k - 1L] <- u
  Sroad <- segments_new(p$x[from_of[er]], p$y[from_of[er]], p$x[sub$target[er]], p$y[sub$target[er]])
  th <- seq(0, 2 * pi, length.out = 37)
  Sbd <- segments_new(R_FLOOD * cos(th[-37]), R_FLOOD * sin(th[-37]), R_FLOOD * cos(th[-1]), R_FLOOD * sin(th[-1]))
  Sroad <- round(Sroad, 6); Sbd <- round(Sbd, 6)         # same 1e-6 grid as the algorithms use
  S <- rbind(Sroad, Sbd); grp <- c(rep(1L, nrow(Sroad)), rep(2L, nrow(Sbd)))
  sw <- sweep_line_intersections(S, group = grp)
  cross_sw <- integer(0)
  for (k in seq_len(nrow(sw$pairs))) { i <- min(sw$pairs[k, ]); if (i <= nrow(Sroad)) cross_sw <- c(cross_sw, i) }
  cross_sw <- unique(cross_sw); cross_sw <- cross_sw[merge_sort_perm(cross_sw)]
  # --- hulls of the incidents inside the flood zone and of the shelters
  pi_ <- local_xy(inc$lon, inc$lat, lon0, lat0); inside <- which(pi_$x^2 + pi_$y^2 <= (1.5 * R_FLOOD)^2)
  hx <- pi_$x[inside]; hy <- pi_$y[inside]
  gh <- graham_scan(hx, hy)
  scen[[s]] <- list(
    id = s, lon0 = lon0, lat0 = lat0, r_flood = R_FLOOD, graph = sub, orig_id = attr(sub, "orig_id"),
    danger = dang, people = people, shelters = shel, shelter_cap = scap, flows = flows, flow_total = pr$value,
    cross_edges = er[cross_sw],
    hull = data.frame(lon = inc$lon[inside][gh], lat = inc$lat[inside][gh]),
    hull_info = list(points = length(inside), hull_size = length(gh), area_km2 = polygon_area(hx[gh], hy[gh]) / 1e6),
    incidents_inside = inside)
}
save_pc(scen, "evacuation")

# =============================================================================
# 3. RESILIENCE (road study area for Karger; max subarray worst load period)
#    Karger itself runs live in the app (fixed 100 trials, fixed seed).
# =============================================================================
R <- study_instance(g, 120L, 7L)
und <- csr_undirected_edges(R$g)
n_hours <- 14 * 24
load <- hourly_counts(inc$time_min, n_hours)
cap_ <- as.numeric(stats::quantile(load, 0.65))
score <- load - cap_
save_pc(list(study = list(graph = R$g, edges = und, orig_id = attr(R$g, "orig_id")),
             load = load, capacity = cap_, window = max_subarray_dc(score)), "resilience")

# =============================================================================
# 4. PLACEMENT / LOGISTICS INSTANCE (set cover, vertex cover; TSP distance matrix)
#    Exact branch and bound runs under a node cap; the app uses it only when it
#    finished (complete = TRUE) and otherwise falls back to the approximation.
# =============================================================================
P <- study_instance(g, 100L, 9L)
U <- P$g$n
qs <- as.numeric(stats::quantile(P$D[upper.tri(P$D)], c(0.10, 0.18, 0.28)))
sc <- list()
for (k in seq_along(qs)) {
  sets <- coverage_sets(P$D, qs[k])
  sc[[k]] <- list(radius_s = qs[k], max_set = max(lengths(sets)), greedy = set_cover_greedy(sets, U),
                  exact = set_cover_exact(sets, U, max_nodes = 60000))
}
pu <- csr_undirected_edges(P$g)
save_pc(list(graph = P$g, D = P$D, orig_id = attr(P$g, "orig_id"), set_cover = sc, edges = pu,
             vertex_cover = list(approx = vertex_cover_approx(pu$n, pu$eu, pu$ev),
                                 exact = vertex_cover_exact(pu$n, pu$eu, pu$ev, max_nodes = 150000))), "placement")

# =============================================================================
# 5. INTEL (LCS de-duplication scores; KMP search runs live in the app)
# =============================================================================
calls <- syn$calls
sel <- c(seq_len(min(450L, sum(!calls$is_duplicate))), which(calls$is_duplicate & calls$dup_of <= 450L))
sub_calls <- calls[sel, ]
save_pc(list(calls = sub_calls, best = report_best_matches(sub_calls$text, 0.5)), "intel")
cat("precompute complete\n")
