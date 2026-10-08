# =============================================================================
# Module 3 - String matching: Knuth-Morris-Pratt (KMP)
# CityShield use: search incident logs / call transcripts for keywords
#   ("flood", "gas cylinder", an incident id ...).
# Texts and patterns are converted ONCE to integer character codes (utf8ToInt),
# so every comparison below is an integer comparison; `comparisons` counts them,
# which is what the theoretical bounds are about.
# kmp_match returns 1-based start positions of ALL (overlapping) matches.
# =============================================================================

#' Characters of a string as integer codes. O(n).
str_codes <- function(s) utf8ToInt(s)

# ---- KMP --------------------------------------------------------------------------

#' KMP failure (prefix) function: pi[j] = length of the longest PROPER prefix of
#' P[1..j] that is also a suffix of it.
#' Approach: each pi[j] extends pi[j-1] by one if the next characters match,
#'   otherwise fall back along earlier borders. The border length k rises by at
#'   most 1 per character, so total fallbacks <= m (amortised).
#' Time complexity: O(m).  Space complexity: O(m).
kmp_prefix <- function(pattern) {
  m <- length(pattern); pi <- integer(m); k <- 0L
  if (m >= 2L) for (j in 2:m) {
    while (k > 0L && pattern[k + 1L] != pattern[j]) k <- pi[k]
    if (pattern[k + 1L] == pattern[j]) k <- k + 1L
    pi[j] <- k
  }
  pi
}

#' Knuth-Morris-Pratt matching.
#'
#' Problem: find every position where pattern P (length m) occurs in text T
#'   (length n), without ever re-examining a text character.
#' Approach: keep q = number of pattern characters currently matched. On a
#'   mismatch do not move back in the text: shrink q to pi[q] (the longest
#'   border) and retry the same text character. The pattern "knows" which
#'   shifts are impossible because of what it already matched.
#' Time complexity: O(n + m) - q grows at most n times in total, so it can fall
#'   at most n times (amortised argument); preprocessing is O(m).
#' Space complexity: O(m) for the prefix table.
#' @return list(pos, comparisons)
kmp_match <- function(text, pattern) {
  n <- length(text); m <- length(pattern)
  pos <- integer(0); comps <- 0
  if (m == 0L || m > n) return(list(pos = pos, comparisons = 0))
  pi <- kmp_prefix(pattern)
  q <- 0L
  for (i in seq_len(n)) {
    while (q > 0L) {
      comps <- comps + 1
      if (pattern[q + 1L] == text[i]) break
      q <- pi[q]
    }
    if (q == 0L) { comps <- comps + 1; if (pattern[1L] == text[i]) q <- 1L }
    else q <- q + 1L
    if (q == m) { pos <- c(pos, i - m + 1L); q <- pi[q] }
  }
  list(pos = pos, comparisons = comps)
}
