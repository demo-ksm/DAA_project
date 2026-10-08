# =============================================================================
# Module 5 - Convex hull: Graham scan
# CityShield use: the smallest convex region containing all affected incidents,
#                 and the coverage hull of a set of stations.
# Returns the STRICT hull (collinear boundary points are dropped) as point
# indices in counter-clockwise order, starting from the lowest-leftmost point.
# =============================================================================

#' Graham scan convex hull.
#'
#' Problem: smallest convex polygon enclosing a point set.
#' Approach: take the lowest point P (leftmost on ties), sort the others by
#'   polar angle around P (ties: keep only the farthest, because nearer
#'   collinear points cannot be vertices), then sweep through them keeping a
#'   stack of hull vertices. A new point q is pushed only after popping every
#'   top vertex that would make a non-left turn (cross <= 0) with q.
#'   Each point is pushed once and popped at most once.
#' Time complexity: O(n log n) - dominated by the angular sort (the scan is O(n)).
#' Space complexity: O(n).
#'
#' @param x,y coordinate vectors
#' @return integer vector of hull indices (CCW)
graham_scan <- function(x, y) {
  n <- length(x)
  if (n == 0L) return(integer(0))
  # pivot = lowest y, then lowest x
  p0 <- 1L
  for (i in seq_len(n)) if (y[i] < y[p0] || (y[i] == y[p0] && x[i] < x[p0])) p0 <- i
  if (n == 1L) return(p0)
  # candidates = all points except the pivot
  idx <- integer(n - 1L); m <- 0L
  for (i in seq_len(n)) if (i != p0) { m <- m + 1L; idx[m] <- i }
  ang <- numeric(m); dst <- numeric(m)
  for (k in seq_len(m)) {
    ang[k] <- atan2(y[idx[k]] - y[p0], x[idx[k]] - x[p0])
    dst[k] <- dist2(x[idx[k]], y[idx[k]], x[p0], y[p0])
  }
  # stable lexicographic sort by (angle, distance): sort by distance first, then angle
  o1 <- merge_sort_perm(dst)
  ang1 <- numeric(m); idx1 <- integer(m); dst1 <- numeric(m)
  for (k in seq_len(m)) { ang1[k] <- ang[o1[k]]; idx1[k] <- idx[o1[k]]; dst1[k] <- dst[o1[k]] }
  o2 <- merge_sort_perm(ang1)
  # keep only the farthest point of each ray from the pivot
  cand <- integer(m); nc <- 0L
  for (k in seq_len(m)) {
    i <- idx1[o2[k]]
    if (nc > 0L && orientation(x[p0], y[p0], x[cand[nc]], y[cand[nc]], x[i], y[i]) == 0L) {
      cand[nc] <- i                      # same ray, later in order = farther
    } else {
      nc <- nc + 1L; cand[nc] <- i
    }
  }
  if (nc == 0L) return(p0)                # all points identical to the pivot
  hull <- integer(nc + 1L); top <- 1L; hull[1L] <- p0
  for (k in seq_len(nc)) {
    i <- cand[k]
    while (top >= 2L &&
           cross(x[hull[top - 1L]], y[hull[top - 1L]], x[hull[top]], y[hull[top]], x[i], y[i]) <= 1e-12) {
      top <- top - 1L                     # not a strict left turn: pop
    }
    top <- top + 1L; hull[top] <- i
  }
  hull[seq_len(top)]
}

#' Polygon area by the shoelace formula (vertices in order). O(h).
polygon_area <- function(x, y) {
  h <- length(x); a <- 0
  if (h < 3L) return(0)
  j <- h
  for (i in seq_len(h)) { a <- a + (x[j] + x[i]) * (y[j] - y[i]); j <- i }
  abs(a) / 2
}

#' Is point (px, py) inside (or on) a CONVEX CCW polygon? O(h).
#' Approach: it must be on the left of every directed edge.
point_in_convex_polygon <- function(px, py, hx, hy) {
  h <- length(hx)
  if (h == 0L) return(FALSE)
  if (h == 1L) return(px == hx[1] && py == hy[1])
  for (i in seq_len(h)) {
    j <- if (i == h) 1L else i + 1L
    if (cross(hx[i], hy[i], hx[j], hy[j], px, py) < -1e-12) return(FALSE)
  }
  TRUE
}
