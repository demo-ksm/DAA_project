# Module 3: Knuth-Morris-Pratt string matching.

# Reference: all (overlapping) occurrences by comparing every substring (test-side only).
ref_matches <- function(text, pattern) {
  n <- nchar(text); m <- nchar(pattern)
  if (m == 0 || m > n) return(integer(0))
  starts <- seq_len(n - m + 1)
  starts[substring(text, starts, starts + m - 1) == pattern]
}
rand_text <- function(n, alphabet, seed) { set.seed(seed); paste(sample(alphabet, n, TRUE), collapse = "") }

test_that("KMP finds all overlapping occurrences (random binary / DNA / English-like text)", {
  for (seed in 1:12) {
    alpha <- list(c("a", "b"), c("a", "c", "g", "t"), c(letters[1:6], " "))[[seed %% 3 + 1]]
    text <- rand_text(sample(30:300, 1), alpha, seed)
    pat <- substr(text, sample(1:10, 1), sample(11:14, 1))            # a pattern that really occurs
    for (p in c(pat, "zz", substr(text, 1, 1))) {
      expect_equal(kmp_match(str_codes(text), str_codes(p))$pos, ref_matches(text, p), info = paste(seed, p))
    }
  }
})

test_that("KMP edge cases: empty pattern, pattern longer than text, equal, single char, repeats", {
  f <- kmp_match
  expect_equal(f(str_codes("abc"), integer(0))$pos, integer(0))
  expect_equal(f(str_codes("ab"), str_codes("abc"))$pos, integer(0))
  expect_equal(f(str_codes("abc"), str_codes("abc"))$pos, 1L)
  expect_equal(f(str_codes("a"), str_codes("a"))$pos, 1L)
  expect_equal(f(str_codes(""), str_codes("a"))$pos, integer(0))
  expect_equal(f(str_codes("aaaaa"), str_codes("aa"))$pos, 1:4)                  # overlaps
  expect_equal(f(str_codes("abababab"), str_codes("abab"))$pos, c(1L, 3L, 5L))
})

test_that("KMP prefix function on classic patterns", {
  expect_equal(kmp_prefix(str_codes("ababaca")), c(0L, 0L, 1L, 2L, 3L, 0L, 1L))
  expect_equal(kmp_prefix(str_codes("aaaa")), c(0L, 1L, 2L, 3L))
  expect_equal(kmp_prefix(str_codes("abcd")), integer(4))
  expect_equal(kmp_prefix(integer(0)), integer(0))
})

test_that("KMP comparison count is linear even on the worst case for naive matching", {
  n <- 2000; m <- 100
  text <- str_codes(strrep("a", n)); pat <- str_codes(paste0(strrep("a", m - 1), "b"))
  km <- kmp_match(text, pat)
  expect_equal(length(km$pos), 0L)
  expect_lte(km$comparisons, 2 * n)                                  # <= 2n (amortised)
  rnd <- str_codes(rand_text(n, letters, 1))
  expect_lte(kmp_match(rnd, str_codes("abc"))$comparisons, 2 * n)
})

test_that("KMP finds keywords in the synthetic call log", {
  skip_if_not(file.exists(data_path("synthetic", "call_logs.rds")))
  cl <- readRDS(data_path("synthetic", "call_logs.rds"))
  txt <- paste(cl$text, collapse = "\n")
  for (kw in c("flood", "fire", "gas cylinder"))
    expect_equal(kmp_match(str_codes(txt), str_codes(kw))$pos, ref_matches(txt, kw), info = kw)
})
