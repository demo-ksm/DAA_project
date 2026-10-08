# Run the whole test suite and store a per-file summary for the report.
#   Rscript scripts/07_test_summary.R
library(testthat)
res <- as.data.frame(test_dir("tests/testthat", reporter = "silent", stop_on_failure = FALSE))
tsum <- do.call(rbind, lapply(split(res, res$file), function(d)
  data.frame(file = d$file[1], tests = nrow(d), expectations = sum(d$nb), failed = sum(d$failed),
             skipped = sum(d$skipped), errors = sum(d$error), warnings = sum(d$warning), seconds = round(sum(d$real), 1))))
rownames(tsum) <- NULL
dir.create("benchmarks/results", showWarnings = FALSE)
saveRDS(tsum, "benchmarks/results/test_summary.rds")
print(tsum)
cat(sprintf("\nTOTAL: %d test blocks, %d expectations, %d failed, %d skipped, %d errors\n",
            sum(tsum$tests), sum(tsum$expectations), sum(tsum$failed), sum(tsum$skipped), sum(tsum$errors)))
