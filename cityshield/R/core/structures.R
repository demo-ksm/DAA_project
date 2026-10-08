# =============================================================================
# core/structures.R  -  Hand-written data structures (no packages)
#
# All structures live in ENVIRONMENTS so push/pop mutate in place (R would
# otherwise copy the whole vector on every modification -> O(n) per operation).
#
#   Queue  : circular buffer, preallocated vector + head/tail indices
#   Stack  : preallocated vector + top index (doubling when full)
#   MinHeap: binary heap in an array; parent(i) = i %/% 2, children 2i, 2i+1
# =============================================================================

# ---- Queue (circular buffer) ------------------------------------------------

#' Fixed-capacity FIFO queue.
#' Problem: O(1) enqueue/dequeue for BFS, FIFO branch & bound.
#' Approach: preallocated integer vector; `head` = index of the front element,
#'   `tail` = index where the next element is written; both wrap modulo capacity.
#' Time: push/pop/empty O(1).  Space: O(capacity).
queue_new <- function(capacity) {
  q <- new.env(parent = emptyenv())
  q$buf <- integer(capacity)
  q$cap <- as.integer(capacity)
  q$head <- 1L
  q$tail <- 1L
  q$size <- 0L
  q
}

queue_empty <- function(q) q$size == 0L
queue_size  <- function(q) q$size

queue_push <- function(q, x) {
  if (q$size == q$cap) stop("queue overflow (capacity ", q$cap, ")")
  q$buf[q$tail] <- x
  q$tail <- if (q$tail == q$cap) 1L else q$tail + 1L
  q$size <- q$size + 1L
  invisible(q)
}

queue_pop <- function(q) {
  if (q$size == 0L) stop("queue underflow")
  x <- q$buf[q$head]
  q$head <- if (q$head == q$cap) 1L else q$head + 1L
  q$size <- q$size - 1L
  x
}

# Same queue but storing doubles (needed when payload is not an integer id).
queue_new_num <- function(capacity) {
  q <- queue_new(1L)
  q$buf <- numeric(capacity)
  q$cap <- as.integer(capacity)
  q
}

# ---- Stack (explicit, for iterative DFS / backtracking) ---------------------

#' Growable LIFO stack of integers.
#' Problem: replace recursion (R's default stack depth ~ 5000 frames is far less
#'   than a road-network DFS path) with an explicit stack.
#' Approach: preallocated vector + `top` index; capacity doubles when full
#'   (amortised O(1) push).
#' Time: push (amortised) / pop / empty O(1).  Space: O(n).
stack_new <- function(capacity = 64L) {
  s <- new.env(parent = emptyenv())
  s$buf <- integer(capacity)
  s$top <- 0L
  s
}

stack_empty <- function(s) s$top == 0L
stack_size  <- function(s) s$top

stack_push <- function(s, x) {
  if (s$top == length(s$buf)) length(s$buf) <- 2L * length(s$buf)
  s$top <- s$top + 1L
  s$buf[s$top] <- x
  invisible(s)
}

stack_pop <- function(s) {
  if (s$top == 0L) stop("stack underflow")
  x <- s$buf[s$top]
  s$top <- s$top - 1L
  x
}

stack_peek <- function(s) {
  if (s$top == 0L) stop("stack empty")
  s$buf[s$top]
}

# ---- Binary min-heap --------------------------------------------------------

#' Binary min-heap keyed by a double, carrying an integer payload.
#' Problem: priority queue for Dijkstra, best-first branch & bound.
#' Approach: complete binary tree stored in arrays `key[]` / `val[]`
#'   (1-based: parent(i) = i %/% 2, children 2i and 2i+1).
#'   push = append + sift-up; pop = move last to root + sift-down.
#'   Capacity doubles when full.
#' Time: push O(log n), pop O(log n), peek O(1).  Space: O(n).
#' (For a MAX-heap push the negated key.)
heap_new <- function(capacity = 64L) {
  h <- new.env(parent = emptyenv())
  h$key <- numeric(capacity)
  h$val <- integer(capacity)
  h$size <- 0L
  h
}

heap_empty <- function(h) h$size == 0L
heap_size  <- function(h) h$size

heap_push <- function(h, key, val) {
  if (h$size == length(h$key)) {
    newcap <- 2L * length(h$key)
    length(h$key) <- newcap
    length(h$val) <- newcap
  }
  i <- h$size + 1L
  h$size <- i
  # sift-up: move the hole upward while the parent is larger
  while (i > 1L) {
    p <- i %/% 2L
    if (h$key[p] <= key) break
    h$key[i] <- h$key[p]
    h$val[i] <- h$val[p]
    i <- p
  }
  h$key[i] <- key
  h$val[i] <- val
  invisible(h)
}

heap_peek_key <- function(h) {
  if (h$size == 0L) stop("heap empty")
  h$key[1L]
}

#' Remove the minimum; returns list(key=, val=).
heap_pop <- function(h) {
  if (h$size == 0L) stop("heap underflow")
  top_key <- h$key[1L]
  top_val <- h$val[1L]
  n <- h$size
  last_key <- h$key[n]
  last_val <- h$val[n]
  n <- n - 1L
  h$size <- n
  if (n > 0L) {
    # sift-down the former last element from the root
    i <- 1L
    repeat {
      c <- 2L * i
      if (c > n) break
      if (c < n && h$key[c + 1L] < h$key[c]) c <- c + 1L   # smaller child
      if (h$key[c] >= last_key) break
      h$key[i] <- h$key[c]
      h$val[i] <- h$val[c]
      i <- c
    }
    h$key[i] <- last_key
    h$val[i] <- last_val
  }
  list(key = top_key, val = top_val)
}

# ---- Binary min-heap with a LEXICOGRAPHIC two-part key ----------------------

#' Min-heap ordered by (k1, then k2). Used as the sweep-line event queue, where
#' events are ordered by x and ties are broken by y.
#' Time: push/pop O(log n).  Space: O(n).
heap2_new <- function(capacity = 64L) {
  h <- new.env(parent = emptyenv())
  h$k1 <- numeric(capacity); h$k2 <- numeric(capacity)
  h$val <- integer(capacity); h$size <- 0L
  h
}
heap2_empty <- function(h) h$size == 0L
heap2_size  <- function(h) h$size

# TRUE when (a1,a2) < (b1,b2) lexicographically
.lex_less <- function(a1, a2, b1, b2) a1 < b1 || (a1 == b1 && a2 < b2)

heap2_push <- function(h, k1, k2, val) {
  if (h$size == length(h$k1)) {
    nc <- 2L * length(h$k1)
    length(h$k1) <- nc; length(h$k2) <- nc; length(h$val) <- nc
  }
  i <- h$size + 1L
  h$size <- i
  while (i > 1L) {
    p <- i %/% 2L
    if (!.lex_less(k1, k2, h$k1[p], h$k2[p])) break
    h$k1[i] <- h$k1[p]; h$k2[i] <- h$k2[p]; h$val[i] <- h$val[p]
    i <- p
  }
  h$k1[i] <- k1; h$k2[i] <- k2; h$val[i] <- val
  invisible(h)
}

heap2_pop <- function(h) {
  if (h$size == 0L) stop("heap underflow")
  top <- list(k1 = h$k1[1L], k2 = h$k2[1L], val = h$val[1L])
  n <- h$size
  l1 <- h$k1[n]; l2 <- h$k2[n]; lv <- h$val[n]
  n <- n - 1L
  h$size <- n
  if (n > 0L) {
    i <- 1L
    repeat {
      c <- 2L * i
      if (c > n) break
      if (c < n && .lex_less(h$k1[c + 1L], h$k2[c + 1L], h$k1[c], h$k2[c])) c <- c + 1L
      if (!.lex_less(h$k1[c], h$k2[c], l1, l2)) break
      h$k1[i] <- h$k1[c]; h$k2[i] <- h$k2[c]; h$val[i] <- h$val[c]
      i <- c
    }
    h$k1[i] <- l1; h$k2[i] <- l2; h$val[i] <- lv
  }
  top
}
