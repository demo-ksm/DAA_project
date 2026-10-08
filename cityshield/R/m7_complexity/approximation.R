# =============================================================================
# Module 7 - Approximation algorithms vs exact solutions (NP-hard problems)
#   Vertex cover  -> SENSOR placement on road junctions (every road link is
#                    watched by a sensor at one of its two ends)
#   Set cover     -> fire/ambulance STATION placement (every demand point is
#                    within reach of a chosen station)
#   TSP 2-approx  -> INSPECTION route through k sites
# Instances are small (50-150 node road subgraphs) because the exact solvers are
# exponential.
# =============================================================================

# ---- Vertex cover ----------------------------------------------------------------

#' Is `cover` a vertex cover of the edge list? O(m).
is_vertex_cover <- function(eu, ev, cover, n) {
  inc <- logical(n); for (v in cover) inc[v] <- TRUE
  for (e in seq_along(eu)) if (!inc[eu[e]] && !inc[ev[e]]) return(FALSE)
  TRUE
}

#' 2-APPROXIMATION for vertex cover (maximal matching).
#' Problem: smallest set of vertices touching every edge (NP-hard).
#' Approach: scan edges; whenever an edge has both endpoints uncovered, take BOTH
#'   endpoints. The chosen edges form a maximal matching M and the cover has
#'   size 2|M|.
#' Guarantee: any cover must contain at least one endpoint of each matched edge
#'   (the matched edges are disjoint), so OPT >= |M| and |cover| = 2|M| <= 2 OPT.
#' Time complexity: O(n + m).  Space complexity: O(n).
vertex_cover_approx <- function(n, eu, ev) {
  inc <- logical(n); size <- 0L; matching <- 0L
  for (e in seq_along(eu)) {
    if (!inc[eu[e]] && !inc[ev[e]]) {
      inc[eu[e]] <- TRUE; inc[ev[e]] <- TRUE; size <- size + 2L; matching <- matching + 1L
    }
  }
  cover <- integer(size); k <- 0L
  for (v in seq_len(n)) if (inc[v]) { k <- k + 1L; cover[k] <- v }
  list(cover = cover, size = size, matching = matching)
}

#' EXACT minimum vertex cover by branch and bound (explicit stack).
#' Approach: take an uncovered edge (u,v): any cover contains u or v, so branch
#'   on both. Reductions: a vertex of degree 1 forces its neighbour; bound:
#'   |chosen| + (uncovered edges / max degree) must beat the incumbent. The
#'   incumbent starts from the 2-approximation.
#' Time complexity: O(2^k) for cover size k (O(1.4^n)-ish with the rules).
#' Space complexity: O(n) per frame, depth <= cover size.
#' @param max_nodes give up (complete = FALSE) after this many expansions
#' @return list(size, cover, nodes, complete)
vertex_cover_exact <- function(n, eu, ev, max_nodes = Inf) {
  m <- length(eu)
  ap <- vertex_cover_approx(n, eu, ev)
  best <- ap$size; best_cover <- ap$cover
  if (m == 0L) return(list(size = 0L, cover = integer(0), nodes = 0L, complete = TRUE))
  # frame: chosen mask (logical n)
  stack <- list(logical(n)); top <- 1L; nodes <- 0L; complete <- TRUE
  while (top > 0L) {
    inc <- stack[[top]]; top <- top - 1L
    nodes <- nodes + 1L
    if (nodes > max_nodes) { complete <- FALSE; break }
    csize <- 0L; for (v in seq_len(n)) if (inc[v]) csize <- csize + 1L
    if (csize >= best) next
    # degrees among uncovered edges and one uncovered edge
    deg <- integer(n); unc <- 0L; fe <- 0L
    for (e in seq_len(m)) if (!inc[eu[e]] && !inc[ev[e]]) {
      unc <- unc + 1L; deg[eu[e]] <- deg[eu[e]] + 1L; deg[ev[e]] <- deg[ev[e]] + 1L
      if (fe == 0L) fe <- e
    }
    if (unc == 0L) {                                        # a cover: improve incumbent
      best <- csize; best_cover <- integer(0); for (v in seq_len(n)) if (inc[v]) best_cover <- c(best_cover, v)
      next
    }
    maxd <- 0L; for (v in seq_len(n)) if (deg[v] > maxd) maxd <- deg[v]
    if (csize + ceiling(unc / maxd) >= best) next            # bound
    # degree-1 rule: a leaf's neighbour must be in some optimal cover
    forced <- 0L
    for (e in seq_len(m)) if (!inc[eu[e]] && !inc[ev[e]]) {
      if (deg[eu[e]] == 1L) { forced <- ev[e]; break }
      if (deg[ev[e]] == 1L) { forced <- eu[e]; break }
    }
    if (forced != 0L) { i2 <- inc; i2[forced] <- TRUE; top <- top + 1L; stack[[top]] <- i2; next }
    # branch on the endpoint of highest degree of the first uncovered edge
    u <- eu[fe]; v <- ev[fe]
    i1 <- inc; i1[u] <- TRUE; i2 <- inc; i2[v] <- TRUE
    top <- top + 1L; stack[[top]] <- if (deg[u] >= deg[v]) i2 else i1
    top <- top + 1L; stack[[top]] <- if (deg[u] >= deg[v]) i1 else i2   # explore the high-degree choice first
  }
  list(size = best, cover = best_cover, nodes = nodes, complete = complete)
}

#' REFERENCE: all subsets in increasing size. Time O(2^n m).  n <= 20.
vertex_cover_brute <- function(n, eu, ev) {
  if (length(eu) == 0L) return(0L)
  best <- n
  for (mask in 0:(2^n - 1)) {
    inc <- logical(n); x <- mask; sz <- 0L
    for (v in seq_len(n)) { if (x %% 2 == 1) { inc[v] <- TRUE; sz <- sz + 1L }; x <- x %/% 2 }
    if (sz >= best) next
    ok <- TRUE
    for (e in seq_along(eu)) if (!inc[eu[e]] && !inc[ev[e]]) { ok <- FALSE; break }
    if (ok) best <- sz
  }
  best
}

# ---- Set cover -------------------------------------------------------------------

#' Greedy set cover (ln n approximation).
#'
#' Problem: universe 1..U, a family of subsets; choose the fewest subsets whose
#'   union is the universe (NP-hard).
#' Approach: repeatedly pick the set covering the MOST still-uncovered elements.
#' Guarantee: the i-th element covered is charged price 1/(new elements of the
#'   chosen set); at that moment at least OPT sets cover the rest, so price
#'   <= OPT/(remaining). Summing gives |greedy| <= OPT * H(d) <= OPT (ln d + 1)
#'   where d is the size of the largest set.
#' Time complexity: O(S * sum of set sizes) simple version (S sets).
#' Space complexity: O(U + S).
#' @param sets list of integer vectors; @param U universe size
#' @return list(chosen, size, covered_all)
set_cover_greedy <- function(sets, U) {
  covered <- logical(U); ncov <- 0L; chosen <- integer(0)
  used <- logical(length(sets))
  while (ncov < U) {
    best <- 0L; bgain <- 0L
    for (s in seq_along(sets)) {
      if (used[s]) next
      gain <- 0L; for (e in sets[[s]]) if (!covered[e]) gain <- gain + 1L
      if (gain > bgain) { bgain <- gain; best <- s }
    }
    if (best == 0L) break                                    # remaining elements uncoverable
    used[best] <- TRUE; chosen <- c(chosen, best)
    for (e in sets[[best]]) if (!covered[e]) { covered[e] <- TRUE; ncov <- ncov + 1L }
  }
  list(chosen = chosen, size = length(chosen), covered_all = ncov == U)
}

#' EXACT set cover by branch and bound (explicit stack, no recursion).
#' Approach: pick the uncovered element with the FEWEST candidate sets; one of
#'   those sets must be in the cover, so branch over exactly them. Bound:
#'   |chosen| + ceil(uncovered / largest remaining set) >= incumbent => prune.
#'   The incumbent starts from the greedy solution.
#' Time complexity: O(f^k) worst case (f = max frequency of an element, k = cover
#'   size).  Space complexity: O(U) per frame, depth <= k.
#' @return list(size, chosen, nodes, complete)
set_cover_exact <- function(sets, U, max_nodes = Inf) {
  gr <- set_cover_greedy(sets, U)
  if (!gr$covered_all) return(list(size = NA_integer_, chosen = integer(0), nodes = 0L, complete = TRUE))
  best <- gr$size; best_sel <- gr$chosen
  S <- length(sets)
  # which sets contain element e
  cont <- vector("list", U); for (e in seq_len(U)) cont[[e]] <- integer(0)
  for (s in seq_len(S)) for (e in sets[[s]]) cont[[e]] <- c(cont[[e]], s)
  fcov <- list(logical(U)); fsel <- list(integer(0)); top <- 1L; nodes <- 0L; complete <- TRUE
  while (top > 0L) {
    cov <- fcov[[top]]; sel <- fsel[[top]]; top <- top - 1L
    nodes <- nodes + 1L
    if (nodes > max_nodes) { complete <- FALSE; break }
    if (length(sel) >= best) next
    unc <- 0L; pick <- 0L; pf <- Inf
    for (e in seq_len(U)) if (!cov[e]) {
      unc <- unc + 1L
      if (length(cont[[e]]) < pf) { pf <- length(cont[[e]]); pick <- e }
    }
    if (unc == 0L) { best <- length(sel); best_sel <- sel; next }
    maxgain <- 0L
    for (s in seq_len(S)) { g <- 0L; for (e in sets[[s]]) if (!cov[e]) g <- g + 1L; if (g > maxgain) maxgain <- g }
    if (length(sel) + ceiling(unc / maxgain) >= best) next   # bound
    for (s in cont[[pick]]) {
      c2 <- cov; for (e in sets[[s]]) c2[e] <- TRUE
      top <- top + 1L; fcov[[top]] <- c2; fsel[[top]] <- c(sel, s)
    }
  }
  list(size = best, chosen = best_sel, nodes = nodes, complete = complete)
}

#' REFERENCE: try all subsets of sets in increasing size. S <= 18.
set_cover_brute <- function(sets, U) {
  S <- length(sets); best <- Inf
  for (mask in 0:(2^S - 1)) {
    sz <- 0L; x <- mask; cov <- logical(U)
    for (s in seq_len(S)) { if (x %% 2 == 1) { sz <- sz + 1L; for (e in sets[[s]]) cov[e] <- TRUE }; x <- x %/% 2 }
    if (sz < best) { all <- TRUE; for (e in seq_len(U)) if (!cov[e]) { all <- FALSE; break }; if (all) best <- sz }
  }
  best
}

#' Coverage sets from a distance matrix: station at node j covers node i when
#' D[j, i] <= radius. O(S * U).
coverage_sets <- function(D, radius) {
  n <- nrow(D); sets <- vector("list", n)
  for (j in seq_len(n)) { s <- integer(0); for (i in seq_len(n)) if (D[j, i] <= radius) s <- c(s, i); sets[[j]] <- s }
  sets
}

#' H(d) = 1 + 1/2 + ... + 1/d, the greedy set-cover approximation factor. O(d).
harmonic_bound <- function(d) { s <- 0; for (i in seq_len(d)) s <- s + 1 / i; s }

# ---- TSP 2-approximation ----------------------------------------------------------

#' Prim's minimum spanning tree on a complete graph with distance matrix D, with
#' a binary heap. Time O(n^2 log n) on the complete graph.  Space O(n).
#' @return list(parent, weight = total MST weight)
mst_prim <- function(D) {
  n <- nrow(D)
  key <- rep(Inf, n); parent <- integer(n); done <- logical(n)
  key[1L] <- 0; h <- heap_new(64L); heap_push(h, 0, 1L); total <- 0
  while (!heap_empty(h)) {
    top <- heap_pop(h); u <- top$val
    if (done[u]) next
    done[u] <- TRUE; total <- total + key[u]
    for (v in seq_len(n)) if (!done[v] && D[u, v] < key[v]) { key[v] <- D[u, v]; parent[v] <- u; heap_push(h, D[u, v], v) }
  }
  list(parent = parent, weight = total)
}

#' TSP 2-approximation for METRIC instances.
#'
#' Problem: short closed tour through all sites (distances obey the triangle
#'   inequality - true for shortest-path distances on the road network).
#' Approach: (1) build an MST; (2) walk it in DFS PREORDER (iterative, explicit
#'   stack) visiting each site once ("shortcutting" repeated visits); (3) close
#'   the tour.
#' Guarantee: MST weight <= OPT (delete one edge of the optimal tour => a
#'   spanning tree). A full walk of the MST costs 2 MST; shortcutting never
#'   increases cost by the triangle inequality, so tour <= 2 MST <= 2 OPT.
#' Time complexity: O(n^2 log n) (MST dominates).  Space complexity: O(n).
#' @return list(tour, cost, mst_weight)
tsp_approx <- function(D) {
  n <- nrow(D)
  if (n == 1L) return(list(tour = 1L, cost = 0, mst_weight = 0))
  mst <- mst_prim(D)
  # children lists
  kids <- vector("list", n); for (v in seq_len(n)) kids[[v]] <- integer(0)
  for (v in 2:n) kids[[mst$parent[v]]] <- c(kids[[mst$parent[v]]], v)
  tour <- integer(n); k <- 0L
  st <- stack_new(n); stack_push(st, 1L)
  while (!stack_empty(st)) {
    u <- stack_pop(st); k <- k + 1L; tour[k] <- u
    ch <- kids[[u]]
    if (length(ch) > 0L) for (i in length(ch):1L) stack_push(st, ch[i])   # preorder
  }
  list(tour = tour, cost = tour_cost(D, tour), mst_weight = mst$weight)
}
