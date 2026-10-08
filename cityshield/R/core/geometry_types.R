# =============================================================================
# core/geometry_types.R  -  Plain data types for Module 5 (geometry)
#
# Points and segments are kept as simple numeric vectors / data frames, so the
# geometry algorithms work on raw coordinates and `sf` is only needed for
# loading data and validating results in tests.
#
#   point   : c(x = , y = )
#   segment : c(x1 = , y1 = , x2 = , y2 = )
#   segment set : numeric matrix, one row per segment, columns x1 y1 x2 y2
# =============================================================================

point <- function(x, y) c(x = x, y = y)

segment <- function(x1, y1, x2, y2) c(x1 = x1, y1 = y1, x2 = x2, y2 = y2)

#' Build a segment matrix from four coordinate vectors.
segments_new <- function(x1, y1, x2, y2) {
  m <- cbind(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
  storage.mode(m) <- "double"
  m
}

#' Cross product (b - a) x (c - a).
#' > 0 : c is LEFT of a->b (counter-clockwise turn)
#' < 0 : c is RIGHT of a->b (clockwise turn)
#' = 0 : a, b, c collinear
#' This single primitive underlies every geometric algorithm in Module 5.
#' Time O(1).
cross <- function(ax, ay, bx, by, cx, cy) {
  (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
}

#' Sign of a number as -1L / 0L / 1L, with an absolute tolerance so floating
#' point noise on (nearly) collinear points is treated as zero.
sgn <- function(v, eps = 1e-12) {
  if (v > eps) 1L else if (v < -eps) -1L else 0L
}

#' Squared Euclidean distance (avoids the sqrt when only comparing). O(1).
dist2 <- function(ax, ay, bx, by) (ax - bx)^2 + (ay - by)^2
