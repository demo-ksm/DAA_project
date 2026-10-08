# =============================================================================
# Module 2 - Dynamic programming: Longest Common Subsequence
# CityShield use: match / de-duplicate incident reports. Two call-log lines
#   describing the same incident share a long common subsequence of WORDS even
#   when words were dropped, re-ordered slightly or changed case.
# =============================================================================

#' LCS length table and one LCS.
#'
#' Problem: longest sequence that is a subsequence of both a and b.
#' Approach (DP): L[i,j] = LCS length of a[1..i], b[1..j].
#'   L[i,j] = L[i-1,j-1] + 1 if a[i] == b[j], else max(L[i-1,j], L[i,j-1]).
#'   Optimal substructure: an LCS of the prefixes ends either with the matching
#'   last symbol or drops one of the last symbols. Trace back from L[m,n].
#' Time complexity: Theta(m n).  Space complexity: Theta(m n) for the table.
#' @param a,b vectors of comparable items (characters, words, integers)
#' @return list(length, subsequence (items), table)
lcs <- function(a, b) {
  m <- length(a); n <- length(b)
  L <- matrix(0L, m + 1L, n + 1L)
  for (i in seq_len(m)) for (j in seq_len(n)) {
    if (a[i] == b[j]) L[i + 1L, j + 1L] <- L[i, j] + 1L
    else L[i + 1L, j + 1L] <- if (L[i, j + 1L] >= L[i + 1L, j]) L[i, j + 1L] else L[i + 1L, j]
  }
  len <- L[m + 1L, n + 1L]
  sub <- a[0]; k <- len
  i <- m; j <- n
  if (len > 0L) sub <- a[seq_len(len)]                # right type, filled below
  while (i > 0L && j > 0L) {
    if (a[i] == b[j]) { sub[k] <- a[i]; k <- k - 1L; i <- i - 1L; j <- j - 1L }
    else if (L[i, j + 1L] >= L[i + 1L, j]) i <- i - 1L
    else j <- j - 1L
  }
  list(length = len, subsequence = sub, table = L)
}

#' LCS length only, in O(min(m,n)) memory (two rolling rows).
#' Time Theta(m n).  Space Theta(min(m, n)).
lcs_length <- function(a, b) {
  if (length(b) > length(a)) { t <- a; a <- b; b <- t }
  n <- length(b)
  prev <- integer(n + 1L); cur <- integer(n + 1L)
  for (i in seq_along(a)) {
    for (j in seq_len(n)) {
      if (a[i] == b[j]) cur[j + 1L] <- prev[j] + 1L
      else cur[j + 1L] <- if (prev[j + 1L] >= cur[j]) prev[j + 1L] else cur[j]
    }
    t <- prev; prev <- cur; cur <- t
  }
  prev[n + 1L]
}

#' REFERENCE: enumerate all 2^m subsequences of `a` and test each against `b`.
#' Time O(2^m * (m + n)).  Space O(m). m <= ~16 only.
lcs_brute <- function(a, b) {
  m <- length(a); best <- 0L
  for (mask in 0:(2^m - 1)) {
    idx <- integer(0); x <- mask
    for (i in seq_len(m)) { if (x %% 2 == 1) idx <- c(idx, i); x <- x %/% 2 }
    k <- length(idx)
    if (k <= best) next
    # is a[idx] a subsequence of b ? (greedy scan)
    p <- 1L
    for (j in seq_along(b)) if (p <= k && b[j] == a[idx[p]]) p <- p + 1L
    if (p > k) best <- k
  }
  best
}

#' Lower-case word tokens of a report (the unit LCS works on). O(len).
report_tokens <- function(text) {
  w <- strsplit(tolower(text), "[^a-z0-9]+")[[1]]
  w[nzchar(w)]
}

#' Remove log-id tokens such as "inc0042" (they identify the LOG LINE, not the
#' incident, so they must not influence similarity). O(len).
strip_report_ids <- function(tokens) tokens[!grepl("^inc[0-9]+$", tokens)]

#' Report similarity = LCS length / length of the longer report, in [0, 1].
report_similarity <- function(a, b) {
  ta <- report_tokens(a); tb <- report_tokens(b)
  if (length(ta) == 0L && length(tb) == 0L) return(1)
  lcs_length(ta, tb) / max(length(ta), length(tb))
}

#' Find duplicate reports: report j is a duplicate of the EARLIEST report i < j
#' whose similarity >= threshold. Tokens are cached; compares only reports
#' inside the same coarse bucket (first INC id is ignored: tokens 2..end).
#' Time O(R^2 * L^2) for R reports of L words (bucketing reduces R).
#' @return integer vector dup_of (NA = original)
dedupe_reports <- function(texts, threshold = 0.75, bucket = NULL) {
  R <- length(texts)
  toks <- vector("list", R)
  for (i in seq_len(R)) toks[[i]] <- strip_report_ids(report_tokens(texts[i]))
  dup <- rep(NA_integer_, R)
  for (j in seq_len(R)) {
    for (i in seq_len(j - 1L)) {
      if (!is.na(dup[i])) next
      if (!is.null(bucket) && bucket[i] != bucket[j]) next
      la <- length(toks[[i]]); lb <- length(toks[[j]])
      if (la == 0L || lb == 0L) next
      if (min(la, lb) / max(la, lb) < threshold) next      # length filter (cheap)
      if (lcs_length(toks[[i]], toks[[j]]) / max(la, lb) >= threshold) { dup[j] <- i; break }
    }
  }
  dup
}

#' For every report j, the most similar EARLIER report i < j (by LCS similarity
#' on word tokens) - computed once so that any duplicate threshold can then be
#' applied in O(R). A cheap length filter skips pairs that cannot reach
#' `min_sim`. Time O(R^2 L^2) worst case.
#' @return data.frame(j, best_i, sim)
report_best_matches <- function(texts, min_sim = 0.5) {
  R <- length(texts)
  toks <- vector("list", R)
  for (i in seq_len(R)) toks[[i]] <- strip_report_ids(report_tokens(texts[i]))
  best_i <- integer(R); sim <- numeric(R)
  for (j in seq_len(R)) {
    lb <- length(toks[[j]])
    for (i in seq_len(j - 1L)) {
      la <- length(toks[[i]])
      if (la == 0L || lb == 0L) next
      if (min(la, lb) / max(la, lb) < min_sim) next
      s <- lcs_length(toks[[i]], toks[[j]]) / max(la, lb)
      if (s > sim[j]) { sim[j] <- s; best_i[j] <- i }
    }
  }
  data.frame(j = seq_len(R), best_i = best_i, sim = sim)
}
