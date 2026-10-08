# Module 5 validation: segments, treap, sweep line, Graham scan (vs sf / grDevices).

pair_keys <- function(m) {                       # order-free comparison of pair matrices
  if (nrow(m) == 0L) return(character(0))
  sort(paste(pmin(m[, 1], m[, 2]), pmax(m[, 1], m[, 2])))
}
# Test-side reference for the sweep line: test every pair of segments (Theta(n^2)).
# Coordinates are snapped to the same 1e-6 grid the sweep line uses.
segments_pairs_ref <- function(S) {
  S <- round(S, 6)
  n <- nrow(S)
  out <- matrix(0L, 0L, 2L)
  if (n < 2L) return(out)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (segments_intersect(S[i, ], S[j, ])) out <- rbind(out, c(i, j))
  out
}
rand_segments <- function(n, seed, span = 100, maxlen = 25) {
  set.seed(seed)
  x1 <- runif(n, 0, span); y1 <- runif(n, 0, span)
  ang <- runif(n, 0, 2 * pi); len <- runif(n, 1, maxlen)
  segments_new(x1, y1, x1 + len * cos(ang), y1 + len * sin(ang))
}
sf_pairs <- function(S) {                        # validation with sf
  lines <- sf::st_sfc(lapply(seq_len(nrow(S)), function(i)
    sf::st_linestring(matrix(S[i, ], 2, 2, byrow = TRUE))))
  m <- sf::st_intersects(lines, sparse = FALSE)
  idx <- which(m & upper.tri(m), arr.ind = TRUE)
  cbind(idx[, 1], idx[, 2])
}

test_that("orientation / on_segment / cross basics", {
  expect_equal(orientation(0, 0, 1, 0, 1, 1), 1L)     # left turn
  expect_equal(orientation(0, 0, 1, 0, 1, -1), -1L)   # right turn
  expect_equal(orientation(0, 0, 1, 1, 2, 2), 0L)
  expect_true(on_segment(0, 0, 4, 4, 2, 2))
  expect_false(on_segment(0, 0, 4, 4, 5, 5))
})

test_that("segment properties: length, angle, parallel, point-segment distance", {
  expect_equal(segment_length(0, 0, 3, 4), 5)
  expect_equal(segment_angle(0, 0, 0, 2), pi / 2)
  expect_true(segments_parallel(c(0, 0, 1, 1), c(0, 1, 2, 3)))
  expect_false(segments_parallel(c(0, 0, 1, 1), c(0, 0, 1, 2)))
  set.seed(1)
  for (i in 1:20) {
    a <- runif(4, 0, 10); p <- runif(2, 0, 10)
    ref <- as.numeric(sf::st_distance(sf::st_point(p), sf::st_linestring(matrix(a, 2, 2, byrow = TRUE))))
    expect_equal(point_segment_distance(p[1], p[2], a), ref)
  }
  expect_equal(point_segment_distance(1, 1, c(0, 0, 0, 0)), sqrt(2))       # degenerate segment
})

test_that("segments_intersect handles all special cases", {
  expect_true(segments_intersect(c(0, 0, 4, 4), c(0, 4, 4, 0)))          # proper cross
  expect_true(segments_intersect(c(0, 0, 4, 0), c(2, 0, 2, 3)))          # T-junction
  expect_true(segments_intersect(c(0, 0, 2, 2), c(2, 2, 5, 0)))          # shared endpoint
  expect_true(segments_intersect(c(0, 0, 4, 0), c(2, 0, 6, 0)))          # collinear overlap
  expect_false(segments_intersect(c(0, 0, 1, 0), c(2, 0, 3, 0)))         # collinear disjoint
  expect_false(segments_intersect(c(0, 0, 1, 1), c(0, 1, 1, 2)))         # parallel
  expect_false(segments_intersect(c(0, 0, 1, 1), c(2, 0, 3, -5)))        # lines meet outside
})

test_that("segment_intersection_point agrees with sf and rejects non-crossings", {
  S <- rand_segments(40, 3)
  for (i in 1:39) for (j in (i + 1):40) {
    p <- segment_intersection_point(S[i, ], S[j, ])
    if (segments_intersect(S[i, ], S[j, ])) {
      expect_false(is.null(p))
      ref <- sf::st_intersection(sf::st_linestring(matrix(S[i, ], 2, 2, byrow = TRUE)),
                                 sf::st_linestring(matrix(S[j, ], 2, 2, byrow = TRUE)))
      expect_equal(unname(p), as.numeric(ref), tolerance = 1e-7)
    } else expect_null(p)
  }
  expect_null(segment_intersection_point(c(0, 0, 1, 1), c(0, 1, 1, 2)))   # parallel
})

test_that("the pairwise reference matches sf::st_intersects", {
  for (seed in 1:3) {
    S <- rand_segments(40, seed)
    expect_equal(pair_keys(segments_pairs_ref(S)), pair_keys(sf_pairs(S)))
  }
  expect_equal(nrow(segments_pairs_ref(segments_new(0, 0, 1, 1))), 0L)
})

# ---- treap -------------------------------------------------------------------

test_that("treap keeps sorted order, supports delete/succ/pred/lower_bound", {
  set.seed(5)
  keys <- sample(1000, 300)                         # distinct keys; id = index
  tr <- treap_new(300)
  for (i in 1:300) treap_insert(tr, i, function(o) keys[i] < keys[o])
  expect_equal(keys[treap_inorder(tr)], sort(keys))
  expect_lte(treap_height(tr), 4 * log2(300))      # random BST: height ~ O(log n)
  gone <- sample(300, 100)
  for (i in gone) treap_delete(tr, i)
  rest <- setdiff(1:300, gone)
  expect_equal(keys[treap_inorder(tr)], sort(keys[rest]))
  ord <- treap_inorder(tr)
  for (k in c(1, 50, 199)) {
    expect_equal(treap_succ(tr, ord[k]), ord[k + 1])
    expect_equal(treap_pred(tr, ord[k + 1]), ord[k])
  }
  expect_equal(treap_succ(tr, ord[200]), 0L); expect_equal(treap_pred(tr, ord[1]), 0L)
  expect_equal(treap_max(tr), ord[200])
  lb <- treap_lower_bound(tr, function(id) keys[id] >= 500)
  expect_equal(keys[lb], min(keys[rest][keys[rest] >= 500]))
  expect_equal(treap_lower_bound(tr, function(id) keys[id] >= 5000), 0L)
  expect_error(treap_delete(tr, gone[1]))
})

test_that("treap edge cases: empty, single, sorted insertion stays shallow", {
  tr <- treap_new(5)
  expect_equal(treap_inorder(tr), integer(0)); expect_equal(treap_max(tr), 0L)
  expect_equal(treap_height(tr), 0L)
  treap_insert(tr, 3L, function(o) TRUE)
  expect_equal(treap_inorder(tr), 3L)
  treap_delete(tr, 3L); expect_equal(tr$size, 0L)
  set.seed(2); tr2 <- treap_new(500)                 # worst case for a plain BST
  for (i in 1:500) treap_insert(tr2, i, function(o) i < o)
  expect_equal(treap_inorder(tr2), 1:500)
  expect_lte(treap_height(tr2), 5 * log2(500))
})

# ---- sweep line --------------------------------------------------------------

test_that("sweep line = pairwise reference on random segments", {
  for (seed in 1:6) {
    S <- rand_segments(60, seed)
    sw <- sweep_line_intersections(S)
    expect_equal(pair_keys(sw$pairs), pair_keys(segments_pairs_ref(S)),
                 info = paste("seed", seed))
  }
})

test_that("sweep line: reported points are real intersection points", {
  S <- rand_segments(50, 11)
  sw <- sweep_line_intersections(S)
  for (k in seq_len(nrow(sw$pairs))) {
    p <- segment_intersection_point(S[sw$pairs[k, 1], ], S[sw$pairs[k, 2], ])
    expect_equal(c(sw$x[k], sw$y[k]), unname(p), tolerance = 1e-6)
  }
})

test_that("sweep line degenerate cases: touching, shared endpoints, vertical, concurrent", {
  # T-junction: endpoint of #2 lies on the interior of #1
  S <- segments_new(c(0, 2), c(0, 0), c(4, 2), c(0, 3))
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), "1 2")
  # chain (shared endpoints): 1-2, 2-3 touch; 1-3 do not
  S <- segments_new(c(0, 1, 2), c(0, 1, 0), c(1, 2, 3), c(1, 0, 1))
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), c("1 2", "2 3"))
  # vertical segment crossing two horizontal ones
  S <- segments_new(c(2, 0, 0), c(-1, 0, 1), c(2, 4, 4), c(2, 0, 1))
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), c("1 2", "1 3"))
  # three lines through one point -> all three pairs
  S <- segments_new(c(-1, -1, 0), c(-1, 1, -1), c(1, 1, 0), c(1, -1, 1))
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), c("1 2", "1 3", "2 3"))
  # star shape: all share the origin as an endpoint
  S <- segments_new(rep(0, 4), rep(0, 4), c(1, -1, 0, 0), c(0, 0, 1, -1))
  expect_equal(nrow(sweep_line_intersections(S)$pairs), 6L)
  # disjoint
  S <- segments_new(c(0, 0), c(0, 5), c(1, 1), c(0, 5))
  expect_equal(nrow(sweep_line_intersections(S)$pairs), 0L)
  # empty / single
  expect_equal(nrow(sweep_line_intersections(matrix(0, 0, 4))$pairs), 0L)
  expect_equal(nrow(sweep_line_intersections(segments_new(0, 0, 1, 1))$pairs), 0L)
})

test_that("sweep line grid with many crossings and vertical/horizontal mix", {
  h <- segments_new(rep(0, 6), (1:6) + 0.5, rep(7, 6), (1:6) + 0.5 + seq(0, 0.05, length.out = 6))
  v <- segments_new((1:6) + 0.3, rep(0, 6), (1:6) + 0.3, rep(8, 6))
  S <- rbind(h, v)
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), pair_keys(segments_pairs_ref(S)))
  expect_equal(nrow(sweep_line_intersections(S)$pairs), 36L)
})

test_that("sweep line group filter keeps only cross-group pairs", {
  S <- rand_segments(40, 4)
  grp <- rep(1:2, each = 20)
  got <- pair_keys(sweep_line_intersections(S, group = grp)$pairs)
  all <- segments_pairs_ref(S)
  keep <- all[grp[all[, 1]] != grp[all[, 2]], , drop = FALSE]
  expect_equal(got, pair_keys(keep))
})

test_that("sweep line event count is about 2n + k and status stays small", {
  S <- rand_segments(80, 8, span = 400, maxlen = 20)          # few crossings
  sw <- sweep_line_intersections(S)
  expect_lte(sw$events, 2 * 80 + nrow(sw$pairs))
  expect_lt(sw$max_status, 80)
})

# ---- convex hulls -------------------------------------------------------------

hull_set <- function(idx) sort(idx)
ref_hull <- function(x, y) {                       # grDevices::chull includes no collinear points
  hull_set(grDevices::chull(x, y))
}

test_that("Graham scan agrees with grDevices::chull", {
  for (seed in 1:8) {
    set.seed(seed); n <- sample(c(5, 30, 200), 1)
    x <- runif(n); y <- runif(n)
    g <- graham_scan(x, y)
    expect_equal(hull_set(g), ref_hull(x, y))
    expect_equal(polygon_area(x[g], y[g]),
                 as.numeric(sf::st_area(sf::st_convex_hull(sf::st_multipoint(cbind(x, y))))))
  }
})

test_that("hull output is counter-clockwise and contains all points", {
  set.seed(9); x <- rnorm(150); y <- rnorm(150)
  h <- graham_scan(x, y)
  for (i in seq_along(h)) {
    j <- if (i == length(h)) 1L else i + 1L
    k <- if (j == length(h)) 1L else j + 1L
    expect_gt(cross(x[h[i]], y[h[i]], x[h[j]], y[h[j]], x[h[k]], y[h[k]]), 0)
  }
  for (i in 1:150) expect_true(point_in_convex_polygon(x[i], y[i], x[h], y[h]))
  expect_false(point_in_convex_polygon(100, 100, x[h], y[h]))
})

test_that("hulls: tiny inputs, duplicates, collinear points, grids, circle", {
  for (f in list(graham_scan)) {
    expect_equal(f(numeric(0), numeric(0)), integer(0))
    expect_equal(f(1, 1), 1L)
    expect_equal(hull_set(f(c(0, 1), c(0, 1))), 1:2)
    expect_equal(hull_set(f(c(0, 1, 2, 3), c(0, 1, 2, 3))), c(1L, 4L))      # all collinear
    expect_equal(length(f(c(0, 0, 0, 1, 1, 1), c(0, 0, 0, 0, 0, 0))), 2L)  # duplicates + collinear
    gx <- rep(0:4, 5); gy <- rep(0:4, each = 5)                              # 5x5 grid: strict hull = 4 corners
    expect_equal(length(f(gx, gy)), 4L)
    t <- seq(0, 2 * pi, length.out = 41)[-41]                                # every point on the hull
    expect_equal(length(f(cos(t), sin(t))), 40L)
    tri <- f(c(0, 4, 2, 2), c(0, 0, 3, 1))                                    # interior point dropped
    expect_equal(hull_set(tri), 1:3)
  }
})

test_that("merge_sort_perm sorts stably and matches order()", {
  set.seed(3)
  for (n in c(0, 1, 2, 7, 100, 1000)) {
    k <- sample(20, n, replace = TRUE) + 0.5
    p <- merge_sort_perm(k)
    expect_equal(k[p], sort(k))
    expect_equal(p, order(k))                       # stable == order() default tie-breaking
  }
})

test_that("heap2 orders lexicographically", {
  h <- heap2_new(2)
  pts <- list(c(2, 1), c(1, 5), c(1, 2), c(3, 0), c(1, 2))
  for (i in seq_along(pts)) heap2_push(h, pts[[i]][1], pts[[i]][2], i)
  out <- character(0)
  while (!heap2_empty(h)) { r <- heap2_pop(h); out <- c(out, paste(r$k1, r$k2)) }
  expect_equal(out, c("1 2", "1 2", "1 5", "2 1", "3 0"))
  expect_error(heap2_pop(h))
})

# ---- regression tests: float-noise degeneracies found on real road-style data -----------
test_that("sweep line: a boundary vertex 6e-14 away from a vertical road still touches it", {
  # minimal failing case found while checking road/flood-boundary data:
  # segment 2 ends at x = 6.1e-14 on the vertical segment 1 (x = 0).
  S <- segments_new(c(0, -173.648177666930), c(884.592, 984.807753012208),
                    c(0, 6.12303176911189e-14), c(1105.74, 1000))
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), "1 2")
  expect_equal(pair_keys(segments_pairs_ref(S)), "1 2")
})

test_that("sweep line = pairwise reference on a road grid crossed by a circular boundary (many touching / collinear / vertical cases)", {
  side <- 7; sp <- 200; segs <- NULL
  for (i in 0:(side - 1)) for (j in 0:(side - 1)) {
    if (j < side - 1) segs <- rbind(segs, c(j * sp, i * sp, (j + 1) * sp, i * sp))     # horizontal edges
    if (i < side - 1) segs <- rbind(segs, c(j * sp, i * sp, j * sp, (i + 1) * sp))     # vertical edges
  }
  segs[, c(1, 3)] <- segs[, c(1, 3)] - 3 * sp; segs[, c(2, 4)] <- segs[, c(2, 4)] - 3 * sp   # centre the grid on (0,0)
  th <- seq(0, 2 * pi, length.out = 37)
  circ <- cbind(500 * cos(th[-37]), 500 * sin(th[-37]), 500 * cos(th[-1]), 500 * sin(th[-1]))   # cos(pi/2) ~ 6e-17, -0 values
  S <- rbind(segs, circ); storage.mode(S) <- "double"
  colnames(S) <- c("x1", "y1", "x2", "y2")
  expect_equal(pair_keys(sweep_line_intersections(S)$pairs), pair_keys(segments_pairs_ref(S)))
  grp <- c(rep(1L, nrow(segs)), rep(2L, nrow(circ)))                  # road-vs-boundary only
  got <- pair_keys(sweep_line_intersections(S, group = grp)$pairs)
  all <- segments_pairs_ref(S)
  expect_equal(got, pair_keys(all[grp[all[, 1]] != grp[all[, 2]], , drop = FALSE]))
})

test_that("sweep line never reports a pair twice (negative-zero event keys)", {
  S <- segments_new(c(0, 0), c(-1105.74, -1000), c(0, 173.6), c(-884.592, -984.8))
  r <- sweep_line_intersections(S)
  expect_equal(nrow(r$pairs), 1L)
})
