# =============================================================================
# Module 5 - Sweep line: report ALL segment intersections (Bentley-Ottmann)
# CityShield use: which road segments cross the flood boundary?
# =============================================================================

#' Report all intersecting pairs of segments with a vertical sweep line.
#'
#' Problem: given n segments, list every intersecting pair (and the point).
#' Approach (Bentley-Ottmann, following de Berg et al.):
#'   * A vertical line sweeps left to right. Events are processed in
#'     lexicographic (x, y) order from a min-heap: segment START points,
#'     segment END points, and INTERSECTION points discovered on the way.
#'   * The STATUS = the segments cut by the line, kept in a treap ordered by
#'     y-coordinate on the line. Two segments can only intersect after they
#'     become NEIGHBOURS, so we test only neighbours, whenever the neighbour
#'     relation changes.
#'   * At event p: find all status segments containing p (a contiguous run found
#'     with a lower-bound search); if p is shared by >= 2 segments (starting
#'     here, ending here or passing through) report all their pairs; delete the
#'     ones ending or passing through p, re-insert those starting or passing
#'     through p ordered by their slope (= their order just AFTER p); test the
#'     new neighbour pairs and schedule new intersection events.
#'   Vertical segments are handled like de Berg's horizontal ones (their
#'   position on the line is the event's y; slope = +Inf puts them last).
#'   Robustness: all coordinates are SNAPPED to a 1e-6 grid first, so points that
#'   are equal up to floating-point noise (e.g. a road end-point lying on a
#'   boundary vertex) become exactly equal; event keys are rounded to
#'   the same resolution as the event-dedupe key (the key only; stored coordinates
#'   stay exact). Assumes no two segments
#'   overlap collinearly (documented limitation).
#' Time complexity: O((n + k) log n) for k reported intersections
#'   (2n endpoint events + k crossing events, each O(log n) treap/heap work).
#'   The brute force is Theta(n^2).
#' Space complexity: O(n + k) (events + status).
#'
#' @param S      numeric matrix n x 4 (x1 y1 x2 y2)
#' @param group  optional integer vector; if given, only pairs from DIFFERENT
#'               groups are reported (e.g. road vs flood boundary)
#' @param seed   seed for treap priorities (deterministic runs)
#' @return list(pairs = k x 2 integer matrix (i<j), x, y = intersection points,
#'              events = number of events processed, max_status = max treap size)
sweep_line_intersections <- function(S, group = NULL, seed = 1L) {
  n <- nrow(S)
  EPS <- 1e-7                                 # tolerance << snap grid (1e-6) ... > rounding error of computed points
  S <- round(S, 6)                            # snap input to the grid (see Robustness above)
  set.seed(seed)
  # ---- normalise: start = lexicographically smaller endpoint (x, then y) ----
  ax <- numeric(n); ay <- numeric(n); bx <- numeric(n); by <- numeric(n)
  for (i in seq_len(n)) {
    if (S[i, 1] < S[i, 3] || (S[i, 1] == S[i, 3] && S[i, 2] <= S[i, 4])) {
      ax[i] <- S[i, 1]; ay[i] <- S[i, 2]; bx[i] <- S[i, 3]; by[i] <- S[i, 4]
    } else {
      ax[i] <- S[i, 3]; ay[i] <- S[i, 4]; bx[i] <- S[i, 1]; by[i] <- S[i, 2]
    }
  }
  vert <- logical(n); slope <- numeric(n)
  for (i in seq_len(n)) {
    vert[i] <- abs(bx[i] - ax[i]) <= EPS
    slope[i] <- if (vert[i]) Inf else (by[i] - ay[i]) / (bx[i] - ax[i])
  }
  # y of segment i on the sweep line x = px (py = y of the current event)
  y_at <- function(i, px, py) {
    if (vert[i]) py else ay[i] + slope[i] * (px - ax[i])
  }
  seg_mat <- function(i) c(ax[i], ay[i], bx[i], by[i])

  # ---- event queue: heap ordered by (x, y); one event per distinct point ----
  # All growing arrays live in the environment E: modifying a closure vector with <<-
  # would copy it on every update (O(size) each), spoiling the O((n+k) log n) bound.
  evheap <- heap2_new(max(64L, 4L * n))
  ev_index <- new.env(hash = TRUE, parent = emptyenv())   # point key -> event number
  E <- new.env(parent = emptyenv())
  E$ev_x <- numeric(2L * n + 16L); E$ev_y <- numeric(2L * n + 16L)
  E$ev_U <- vector("list", 2L * n + 16L)                  # segments STARTING at each event
  E$n_ev <- 0L
  E$rp <- integer(16L); E$rpj <- integer(16L); E$rx <- numeric(16L); E$ry <- numeric(16L); E$nres <- 0L
  point_key <- function(x, y) sprintf("%.8f|%.8f", round(x, 8) + 0, round(y, 8) + 0)   # +0 turns -0 into 0
  get_event <- function(x, y) {                           # find or create
    key <- point_key(x, y)    # only the KEY is rounded: stored coordinates stay exact (steep slopes!)
    idx <- ev_index[[key]]
    if (is.null(idx)) {
      E$n_ev <- E$n_ev + 1L
      if (E$n_ev > length(E$ev_x)) {
        cap <- 2L * length(E$ev_x)
        length(E$ev_x) <- cap; length(E$ev_y) <- cap; length(E$ev_U) <- cap
      }
      E$ev_x[E$n_ev] <- x + 0; E$ev_y[E$n_ev] <- y + 0; E$ev_U[[E$n_ev]] <- integer(0)
      assign(key, E$n_ev, envir = ev_index)
      heap2_push(evheap, x, y, E$n_ev)
      idx <- E$n_ev
    }
    idx
  }
  for (i in seq_len(n)) {
    e1 <- get_event(ax[i], ay[i]); E$ev_U[[e1]] <- c(E$ev_U[[e1]], i)
    get_event(bx[i], by[i])
  }

  status <- treap_new(n)
  add_result <- function(i, j, x, y) {
    k <- E$nres + 1L; E$nres <- k
    if (k > length(E$rp)) {
      cap <- 2L * length(E$rp)
      length(E$rp) <- cap; length(E$rpj) <- cap; length(E$rx) <- cap; length(E$ry) <- cap
    }
    E$rp[k] <- min(i, j); E$rpj[k] <- max(i, j); E$rx[k] <- x; E$ry[k] <- y
  }

  # schedule the intersection of two (new) neighbours if it lies after p
  check_pair <- function(s1, s2, px, py) {
    if (s1 == 0L || s2 == 0L) return(invisible())
    q <- segment_intersection_point(seg_mat(s1), seg_mat(s2))
    if (is.null(q)) return(invisible())
    after <- q[1] > px + EPS || (abs(q[1] - px) <= EPS && q[2] > py + EPS)
    if (after) get_event(q[1], q[2])
    invisible()
  }

  n_events <- 0L; max_status <- 0L
  while (!heap2_empty(evheap)) {
    top <- heap2_pop(evheap)
    px <- top$k1; py <- top$k2; idx <- top$val
    n_events <- n_events + 1L
    U <- E$ev_U[[idx]]

    # (1) status segments containing p: contiguous run starting at lower bound
    run <- integer(0)
    cur <- treap_lower_bound(status, function(id) y_at(id, px, py) >= py - EPS)
    while (cur != 0L && y_at(cur, px, py) <= py + EPS) {
      run <- c(run, cur)
      cur <- treap_succ(status, cur)
    }
    Lset <- integer(0); Cset <- integer(0)
    for (s in run) {
      if (abs(bx[s] - px) <= EPS && abs(by[s] - py) <= EPS) Lset <- c(Lset, s)
      else Cset <- c(Cset, s)
    }
    # (2) report all pairs among U u L u C
    all_s <- c(U, Lset, Cset)
    m <- length(all_s)
    if (m >= 2L) {
      for (a in seq_len(m - 1L)) for (b in (a + 1L):m) {
        i <- all_s[a]; j <- all_s[b]
        if (is.null(group) || group[i] != group[j]) add_result(i, j, px, py)
      }
    }
    # (3) remove ending / passing segments; (4) insert starting / passing ones
    for (s in run) treap_delete(status, s)
    for (s in c(U, Cset)) {
      ys <- y_at(s, px, py); sl <- slope[s]
      goes_left <- function(o) {
        yo <- y_at(o, px, py)
        if (ys < yo - EPS) TRUE
        else if (ys > yo + EPS) FALSE
        else sl < slope[o]
      }
      treap_insert(status, s, goes_left)
    }
    if (status$size > max_status) max_status <- status$size
    # (5) new neighbour pairs
    if (length(U) + length(Cset) == 0L) {
      sr <- treap_lower_bound(status, function(id) y_at(id, px, py) >= py - EPS)
      sl_ <- if (sr == 0L) treap_max(status) else treap_pred(status, sr)
      check_pair(sl_, sr, px, py)
    } else {
      lowest <- treap_lower_bound(status, function(id) y_at(id, px, py) >= py - EPS)
      check_pair(treap_pred(status, lowest), lowest, px, py)
      highest <- lowest
      repeat {
        nx <- treap_succ(status, highest)
        if (nx == 0L || y_at(nx, px, py) > py + EPS) break
        highest <- nx
      }
      check_pair(highest, treap_succ(status, highest), px, py)
    }
  }
  k <- E$nres
  list(pairs = cbind(E$rp[seq_len(k)], E$rpj[seq_len(k)]),
       x = E$rx[seq_len(k)], y = E$ry[seq_len(k)],
       events = n_events, max_status = max_status)
}
