test_that("Sumstats are correctly processed for plotting", {
  dat <- data.table::data.table(
    chromosome = c(1, 1, 2, 3),
    base_pair_location = c(1, 2, 3, 4),
    p_value = c(0.1, 0.2, 0.3, 0.4)
  )

  sumstats <- process_sumstats_for_manhattan(dat)

  expect_equal(sumstats$data$chr, c(1, 1, 2, 3))
  expect_equal(sumstats$data$bp_cum, c(3, 4, 8, 13))
  expect_equal(sumstats$data$bp_add[1:4], c(2, 2, 5, 9))
})
