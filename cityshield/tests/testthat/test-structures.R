test_that("queue is FIFO and wraps around the circular buffer", {
  q <- queue_new(3)
  expect_true(queue_empty(q))
  queue_push(q, 1L); queue_push(q, 2L); queue_push(q, 3L)
  expect_error(queue_push(q, 4L), "overflow")
  expect_equal(queue_pop(q), 1L)
  queue_push(q, 4L)                       # reuses the freed slot (wrap-around)
  expect_equal(c(queue_pop(q), queue_pop(q), queue_pop(q)), c(2L, 3L, 4L))
  expect_true(queue_empty(q))
  expect_error(queue_pop(q), "underflow")
})

test_that("queue handles interleaved push/pop over many wraps", {
  q <- queue_new(5)
  expected <- 1L
  for (i in 1:1000) {
    queue_push(q, i)
    if (i %% 2 == 0) { expect_equal(queue_pop(q), expected); expected <- expected + 1L }
    if (queue_size(q) == 5L) { expect_equal(queue_pop(q), expected); expected <- expected + 1L }
  }
})

test_that("stack is LIFO and grows past its initial capacity", {
  s <- stack_new(2)
  expect_true(stack_empty(s))
  for (i in 1:100) stack_push(s, i)
  expect_equal(stack_size(s), 100L)
  expect_equal(stack_peek(s), 100L)
  out <- integer(100)
  for (i in 1:100) out[i] <- stack_pop(s)
  expect_equal(out, 100:1)
  expect_error(stack_pop(s), "underflow")
})

test_that("heap pops keys in non-decreasing order (vs sort() baseline)", {
  set.seed(1)
  for (n in c(0, 1, 2, 10, 500)) {
    keys <- sample.int(1000, n, replace = TRUE) + runif(n)
    h <- heap_new(2)                      # tiny capacity forces growth
    for (i in seq_len(n)) heap_push(h, keys[i], i)
    got <- numeric(n); payload <- integer(n)
    for (i in seq_len(n)) { r <- heap_pop(h); got[i] <- r$key; payload[i] <- r$val }
    expect_equal(got, sort(keys))         # sort() = reference baseline only
    expect_equal(keys[payload], got)      # payload stayed attached to its key
    expect_true(heap_empty(h))
  }
})

test_that("heap handles duplicates, negative keys, and interleaved ops", {
  h <- heap_new()
  heap_push(h, 5, 1L); heap_push(h, 5, 2L); heap_push(h, -3, 3L); heap_push(h, 0, 4L)
  expect_equal(heap_peek_key(h), -3)
  expect_equal(heap_pop(h)$val, 3L)
  heap_push(h, -10, 5L)
  expect_equal(heap_pop(h)$val, 5L)
  expect_equal(heap_pop(h)$key, 0)
  expect_equal(heap_pop(h)$key, 5)
  expect_equal(heap_pop(h)$key, 5)
  expect_error(heap_pop(h), "underflow")
  expect_error(heap_peek_key(h), "empty")
})
