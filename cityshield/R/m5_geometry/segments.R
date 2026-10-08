# =============================================================================
# Module 5 - Line segment properties and segment intersection
# CityShield use: does a road cross the flood boundary? how far is a station
#                 from a road? (all built on the single `cross` primitive)
# =============================================================================

#' Orientation of the ordered triple (p, q, r).
#' Problem: which way do we turn going p -> q -> r?
#' Approach: sign of the cross product (q-p) x (r-p), with an epsilon.
#' Time O(1), space O(1).
#' @return 1 = counter-clockwise (left turn), -1 = clockwise, 0 = collinear
orientation <- function(px, py, qx, qy, rx, ry) {
  sgn(cross(px, py, qx, qy, rx, ry))
}

#' Is point r on segment pq (given that p, q, r are collinear)? O(1).
on_segment <- function(px, py, qx, qy, rx, ry, eps = 1e-12) {
  rx >= min(px, qx) - eps && rx <= max(px, qx) + eps &&
    ry >= min(py, qy) - eps && ry <= max(py, qy) + eps
}

#' Length of a segment. O(1).
segment_length <- function(x1, y1, x2, y2) sqrt(dist2(x1, y1, x2, y2))

#' Direction angle (radians, -pi..pi) of the segment from (x1,y1) to (x2,y2). O(1).
segment_angle <- function(x1, y1, x2, y2) atan2(y2 - y1, x2 - x1)

#' Are two segments parallel? (cross product of their direction vectors = 0). O(1).
segments_parallel <- function(a, b, eps = 1e-12) {
  abs((a[3] - a[1]) * (b[4] - b[2]) - (a[4] - a[2]) * (b[3] - b[1])) <= eps
}

#' Distance from point (px, py) to segment a = c(x1,y1,x2,y2). O(1).
#' Approach: project the point on the supporting line, clamp to the segment.
point_segment_distance <- function(px, py, a) {
  dx <- a[3] - a[1]; dy <- a[4] - a[2]
  len2 <- dx * dx + dy * dy
  if (len2 == 0) return(sqrt(dist2(px, py, a[1], a[2])))
  t <- ((px - a[1]) * dx + (py - a[2]) * dy) / len2
  if (t < 0) t <- 0 else if (t > 1) t <- 1
  sqrt(dist2(px, py, a[1] + t * dx, a[2] + t * dy))
}

#' Do segments a and b (c(x1,y1,x2,y2)) intersect? (CLRS ANY-SEGMENTS test)
#'
#' Problem: decide whether two closed segments share at least one point.
#' Approach: compute the four orientations. If p3,p4 lie on opposite sides of
#'   line p1p2 AND p1,p2 lie on opposite sides of line p3p4 the segments cross
#'   ("general case"). Otherwise they intersect only if some endpoint is
#'   collinear with and lies ON the other segment ("special cases").
#' Time complexity: O(1).  Space complexity: O(1).
segments_intersect <- function(a, b) {
  o1 <- orientation(a[1], a[2], a[3], a[4], b[1], b[2])
  o2 <- orientation(a[1], a[2], a[3], a[4], b[3], b[4])
  o3 <- orientation(b[1], b[2], b[3], b[4], a[1], a[2])
  o4 <- orientation(b[1], b[2], b[3], b[4], a[3], a[4])
  if (o1 != o2 && o3 != o4) return(TRUE)
  if (o1 == 0L && on_segment(a[1], a[2], a[3], a[4], b[1], b[2])) return(TRUE)
  if (o2 == 0L && on_segment(a[1], a[2], a[3], a[4], b[3], b[4])) return(TRUE)
  if (o3 == 0L && on_segment(b[1], b[2], b[3], b[4], a[1], a[2])) return(TRUE)
  if (o4 == 0L && on_segment(b[1], b[2], b[3], b[4], a[3], a[4])) return(TRUE)
  FALSE
}

#' Intersection POINT of two non-parallel segments, or NULL if they do not meet
#' in exactly one point (parallel / collinear / disjoint). O(1).
#' Approach: solve a1 + t*da = b1 + u*db with Cramer's rule; accept 0<=t,u<=1.
segment_intersection_point <- function(a, b, eps = 1e-12) {
  dax <- a[3] - a[1]; day <- a[4] - a[2]
  dbx <- b[3] - b[1]; dby <- b[4] - b[2]
  den <- dax * dby - day * dbx
  if (abs(den) <= eps * (abs(dax) + abs(day) + 1) * (abs(dbx) + abs(dby) + 1)) return(NULL)
  t <- ((b[1] - a[1]) * dby - (b[2] - a[2]) * dbx) / den
  u <- ((b[1] - a[1]) * day - (b[2] - a[2]) * dax) / den
  tol <- 1e-9
  if (t < -tol || t > 1 + tol || u < -tol || u > 1 + tol) return(NULL)
  c(x = a[1] + t * dax, y = a[2] + t * day)
}
